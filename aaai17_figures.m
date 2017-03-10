%% generates figures for the AAAI 2017 paper

% calibrates using the datasets taken on 2/26/16, 7/26/16, and 11/2/16
% loads the dataset taken on 11/1/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

DATADIR = '/Volumes/shared/Projects/Proton Pack/Data';

%% sphere calibration (see go_sphere_again.m)

% vicon data
x = csvload([DATADIR '/20160726/sphere/4/vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

[c, r, e, cs] = sphereFit_ransac(x(:,2:4), 50); % FIXME allowing a huge amount of noise to get a reasonable number of inliers

d = nan([size(x,1) 3]);
for i = find(cs)
    cc = xfconv(x(i,2:7)) \ [c 1]';
    d(i,:) = cc(1:3);
end
d = nanmean(d);

% the product of sphere calibration is d
spherecalib.x = x;
spherecalib.c = c;
spherecalib.r = r;
spherecalib.cs = cs;
spherecalib.d = d;

%% free-space calibration (see go_bias.m)

date = '20160226';
offset = 22.848; % manual (from taps)
material = 'free';
tool = 'stick';
rep = '1';

set = [material rep tool];
prefix = [DATADIR filesep date filesep set filesep];

v = csvload([prefix 'vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

int = csvload([prefix 'teensy.ft.csv'], ...
              [{'Timestamp'}, ...
               arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
           
int = process_mini40(100*ones(size(int,1),1), int, zeros(1,6));

a = round(size(int,1)*1/5);
b = round(size(int,1)*4/5);
int_s = int(a:b,:);

[~,aa] = min(abs(v(:,1)-int(a,1)));
[~,bb] = min(abs(v(:,1)-int(b,1)));
v_s = v(aa:bb,:);

[c,r,e,inliers] = sphereFit_ransac(int_s(:,2:4)); % 94% inliers

[~, start] = min(abs(int_s(:,1) - v(1,1)));

R_m402bod = [ 0       0       1
              1       0       0
              0       1       0 ];

grav = zeros(size(int_s,1), 4);
grav(:,1) = int_s(:,1) - offset;
ideal_grav = grav;
for tm=1:size(grav,1)
    [delta, tv] = min(abs(v(:,1) - grav(tm,1)));
    assert(delta < 0.02);
    Rvw = xfconv(v(tv,5:7));
    grav(tm,2:4) = R_m402bod * int_s(tm,2:4)';
    ideal_grav(tm,2:4) = Rvw\[0 0 -r]';
end

[R, t, err] = rigid_ransac(ideal_grav(inliers,2:4), grav(inliers,2:4), 10); % 100% inliers??

% the product of free-space calibration is R
freecalib.int_s = int_s;
freecalib.grav = grav;
freecalib.ideal_grav = ideal_grav;
freecalib.c = c;
freecalib.r = r;
freecalib.cs = inliers;
freecalib.R = R;
freecalib.t = t;

%% end-effector weighing

date = '20161102';
[~,f] = load_stick([DATADIR filesep date filesep 'weigh/1/']);
[c,r] = sphereFit_ransac(f(:,2:4)); % 98.9% inliers

mass = r/9.81;

%% process calibrations: use d+R to calculate x+y+z for H_vic2bod + H_bal2imu

H_vic2bod = [ 0.9912 -0.0238 -0.1302    0
              0.0164  0.9982 -0.0575   36
              0.1314  0.0549  0.9898 -511.32
              0       0       0         1   ];
H_m402bod = [ 0       0       1       108.99
              1       0       0         0.53
              0       1       0        -2.98
              0       0       0         1   ];
H_bal2imu = [ 1       0       0       254.84
              0       1       0         0
              0       0       1         0
              0       0       0         1   ];

syms x y z real;
R = freecalib.R';
H1 = [R [0 y z]';      0 0 0 1];
H2 = [eye(3) [x 0 0]'; 0 0 0 1];
S = solve(H1 * H2 * [0 0 0 1]' == [d 1]');

H_vic2bod(1:3,1:3) = R;
H_vic2bod(2,4) = S.y;
H_vic2bod(3,4) = S.z;
H_bal2imu(1,4) = S.x;

%% setup for learning

data120 = icra17_load(DATADIR, '20161101', 'stickcam', @(n) n ~= 2);
data120 = icra17_process('bluefox', data120, mass, H_vic2bod, H_m402bod, H_bal2imu);

%% learning

d = data120('vinyl');
start = 17117;
stop = 362378;
t_op = tic;
[vchunks, vp, verr, vX, vY] = online_prediction(d, start, stop, false);
toc(t_op)

% find threshold
stop_point = find(mean(movingmean(abs(diff(vp) ./ vp(2:end,:)), 10), 2) < 1e-3, 1)

%% figures

% COLORS
% red = normal force
% blue = tangential force
% green = tangential speed
% magenta = incremental error
% cyan = overall error
% black = stopping point

% unity-gain first-order low-pass filter
lpfa = .01;
lpfb = [1 lpfa-1];
lpf = @(x) filtfilt(lpfa, lpfb, x);

force = [d.biws(start:stop,1) lpf(d.biws(start:stop,2:end))];
pos   = [d.bvei(start:stop,1) lpf(d.bvei(start:stop,2:end))];
speed = [d.bvei(start:stop,1) [lpf(diff(pos(:,2:end)) / mean(diff(pos(:,1)))); zeros(1,size(pos,2)-1)]];

Ft = sum(sqrt(force(:,2:3).^2),2);
%Ft = zeros(size(force,1), 1);
%for i=1:length(Ft)
%    Ft(i) = dot(force(i,2:3), -speed(i,2:3)/sqrt(sum(speed(i,2:3).^2)));
%end
%Ft(Ft<0) = 0;

% I/O figure
fig_io = figure;
subplot(211);
plot(force(:,1)-force(1,1), force(:,4), 'r', ...
     force(:,1)-force(1,1), Ft, 'b');
xlabel('Time (s)', 'fontsize',16);
ylabel('Force (N)', 'fontsize',16);
set(legend('Normal force', 'Tangential force'), 'fontsize',16);
subplot(212);
plot(speed(:,1)-speed(1,1), sqrt(sum(speed(:,2:3).^2,2)), 'g');
xlabel('Time (s)', 'fontsize',16);
ylabel('Speed (mm/s)', 'fontsize',16);
set(legend('Tangential speed'), 'fontsize',16);

% error figure
fig_error = figure;
subplot(211);
plot(verr(:,1), 'm');
ax = axis;
ax(3) = 0;
axis(ax);
line([stop_point stop_point], ax(3:4), 'color','k');
ylabel('Incremental Error (N)', 'color','k', 'fontsize',16);
subplot(212);
plot(verr(:,2), 'c');
ax = axis;
ax(3) = 0;
axis(ax);
line([stop_point stop_point], ax(3:4), 'color','k');
xlabel('Iteration number', 'fontsize',16);
ylabel('Overall Error (N)', 'fontsize',16);

% coefficient figure
fig_coeff = figure;
subplot(211);
plot(vp(:,1), 'r');
ax = axis;
ax(3) = 0;
axis(ax);
line([stop_point stop_point], ax(3:4), 'color','k');
ylabel('Fn coefficient (N/N)', 'fontsize',16);
subplot(212);
plot(vp(:,2), 'b');
ax = axis;
ax(3) = 0;
axis(ax);
line([stop_point stop_point], ax(3:4), 'color','k');
xlabel('Iteration number', 'fontsize',16);
ylabel('Vt coefficient (Ns/mm)', 'fontsize',16);

% model figure
%% fig_model = figure;
allin = cell2mat(arrayfun(@(s) s.in', vchunks, 'uniformoutput',false))';
allout = reshape([vchunks.out], [],1);
allt = reshape([vchunks.t], [],1) - vchunks(1).t(1);
plot(allt, allout, ...
     allt, allout - allin*vp(stop_point,:)', ...
     allt, allout - allin*vp(end,:)');
xlabel('Time (s)', 'fontsize',16);
ylabel('Tangential force (N)', 'fontsize',16);
set(legend('Actual value', 'Prediction error by model at stopping point', 'Prediction error by model with all data'), 'fontsize',16);
fprintf('RMS Prediction errors: stopping point %.2f, all data %.2f\n', sqrt(mean((allout - allin*vp(stop_point,:)').^2)), sqrt(mean((allout - allin*vp(end,:)').^2)));

%% histogram figure
fig_fhist = figure;
fhistfull = histogram(vX(:,1));
fig_vhist = figure;
vhistfull = histogram(vX(:,2));
fhists = zeros(fhistfull.NumBins, 10);
vhists = zeros(vhistfull.NumBins, 10);
for i=1:10
    frac = 1:i*floor(size(vX,1)/10);
    fhists(:,i) = histcounts(vX(frac,1), fhistfull.BinEdges);
    vhists(:,i) = histcounts(vX(frac,2), vhistfull.BinEdges);
end
figure(fig_fhist);
bar3(fhists(:, [2 5 10]));
view(-90, 0);
ylabel('Normal force (N)', 'fontsize',16);
zlabel('Incidence', 'fontsize',16);
title('Normal force histograms', 'fontsize',16);
set(legend('24 seconds', '60 seconds', '120 seconds', 'location',[.8 .75 0 0]), 'fontsize',16);
figure(fig_vhist);
bar3(vhists(:, [2 5 10]));
view(-90, 0);
ylabel('Tangential speed (mm/s)', 'fontsize',16);
zlabel('Incidence', 'fontsize',16);
title('Tangential speed histograms', 'fontsize',16);
set(legend('24 seconds', '60 seconds', '120 seconds', 'location',[.8 .75 0 0]), 'fontsize',16);

%% save plots
figure(fig_io);
print -dpdf aaai17_io.pdf
figure(fig_error);
print -dpdf aaai17_error.pdf
figure(fig_coeff);
print -dpdf aaai17_coeff.pdf
figure(fig_model);
print -dpdf aaai17_model.pdf
figure(fig_fhist);
print -dpdf aaai17_fhist.pdf
figure(fig_vhist);
print -dpdf aaai17_vhist.pdf
