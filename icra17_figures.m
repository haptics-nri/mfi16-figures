%% generates figures for the ICRA 2017 paper

% calibrates using the datasets taken on 2/26/16, 7/26/16, and 8/11/16
% loads the dataset taken on 8/15/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')
addpath(genpath('geom3d'));

DATADIR = '../data';

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

date = '20160811';
load mini40_calib
[mass, fbias, ~, com, tbias] = weigh({fullfile(DATADIR, date, 'weigh', '1')});
tbias = tbias';
save mini40_calib fbias tbias;

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

% dataset parameters
date = '20160906';
flowtype = 'stickcam';

% end-effector mass comes from calibration above

[data14, materials14] = icra17_load(DATADIR, date, flowtype, @(x) x < 7);
[data38, materials38] = icra17_load(DATADIR, date, flowtype, @(x) x > 5);
assert(all(cellfun(@strcmp, materials14, materials38)));
materials = materials14;

%% load in some manually picked values

% NB to find these values:
% 1. plot(d.v(:,1)-d.v(1,1), d.v(:,2:4), d.int(:,1)-d.v(1,1), d.int(:,2:4))
% 2. d.off = -mean([force spike maxima] - [position spike minima]);
% 3. run the icra17_process steps
% 4. plot(d.iws(:,2:4))
% 5. d.ss = [index just after second spike, index just before third spike];
% 6. plot(d.biws(:,2:4))
% 7. d.bss = [index just after second spike, index just before third spike];

d = data14('abs');
    d.off = 16.0213;
    d.ss =  [28490 125800];
    d.bss = [15400 113000];
data14('abs') = d;
d = data14('glitter');
    d.off = 16.0197;
    d.ss =  [32500 121400];
    d.bss = [20130 108500];
data14('glitter') = d;
d = data14('silk');
    d.off = 16.0192;
    d.ss =  [31560 124600];
    d.bss = [18710 111500];
data14('silk') = d;
d = data14('vinyl');
    d.off = 16.0197;
    d.ss =  [27220 131200];
    d.bss = [14260 118600];
data14('vinyl') = d;
d = data14('wood');
    d.off = 16.0158;
    d.ss =  [27520 128800];
    d.bss = [14210 116000];
data14('wood') = d;

d = data38('abs');
    d.off = 15.9987;
    d.ss =  [24647 140441];
    d.bss = [11935 127889];
data38('abs') = d;
d = data38('glitter');
    d.off = 16.0027;
    d.ss =  [30130 124500];
    d.bss = [17800 111800];
data38('glitter') = d;
d = data38('silk');
    d.off = 16.0005;
    d.ss =  [28960 121400];
    d.bss = [17220 108900];
data38('silk') = d;
d = data38('vinyl');
    d.off = 16.0048;
    d.ss =  [26920 130700];
    d.bss = [14470 118200];
data38('vinyl') = d;
d = data38('wood');
    d.off = 16.0030;
    d.ss =  [29670 124600];
    d.bss = [16840 111800];
data38('wood') = d;

%% process

data14 = icra17_process('both', data14, mass, H_vic2bod, H_m402bod, H_bal2imu);
data38 = icra17_process('both', data38, mass, H_vic2bod, H_m402bod, H_bal2imu);

%% SVM stuff (following Romano & KJK 2014 + Strese & Schuwerk & Steinbach 2015)

[features14, split_idx14, bfeatures14, bsplit_idx14] = icra17_svm(data14, mass);
[features38, split_idx38, bfeatures38, bsplit_idx38] = icra17_svm(data38, mass);

%% grid searches!

all_grid;

%% confusion matrices -- first set gsi to optimal and run the test set

fig_confusion(conf14vp, {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'}, 14, 'Arial', 0, 0.15, true);
print -dpdf icra17_confusion_precision_vicon.pdf;
fig_confusion(conf38bn, {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'}, 14, 'Arial', 0, 0.15, true);
print -dpdf icra17_confusion_precision_bluefox.pdf;

%% accelerometer comparison figure

[v,f,da,dg,~,a,~,dt,~,~,m] = load_stick(fullfile(DATADIR, '20160906', 'stickcam', '6'));
t = cumsum(bitand(dt, 65535))/1e6;
figure
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(211)
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(212)
plot(t, bsxfun(@minus, a(:,5:7), mean(a(:,5:7))))
subplot(211)
axis([9 17 -8 8])
subplot(212)
axis([9 17 -8 8])
xlabel('Time (s)', 'FontSize',16)
ylabel({'Internal accelerometer' 'signal (m/s^2)'}, 'FontSize',16)
subplot(211)
xlabel('Time (s)', 'FontSize',16)
ylabel({'External accelerometer' 'signal (m/s^2)'}, 'FontSize',16)
print -dpdf -r0 icra17_accel_compare.pdf

%% motrak comparison figure

[v,~,~,~,~,~,~,~,~,~,m] = load_stick(fullfile(DATADIR, '20160906', 'stickcam', '6'));
off = [16.028 1814 127]; % 20160906-6
%off = [7.825 817 34]; % 20160712-1
tf = xfconv(v(off(2),2:7)) / xfconv(m(off(3),2:7));
for i=1:size(m,1)
    m(i,2:7) = xfconv(tf * xfconv(m(i,2:7)));
    m(i,[7 6 5]) = rotation3dToEulerAngles(xfconv(m(i,5:7)))*pi/180;
end
for i=1:size(v,1)
    v(i,[7 6 5]) = rotation3dToEulerAngles(xfconv(v(i,5:7)))*pi/180;
end
v(:,5:7) = filter(ones(5,1)/5, 1, v(:,5:7));
subplot(211);
plot(v(:,1) - v(1,1), v(:,2:4), m(:,1) - v(1,1) + off(1), m(:,2:4));
legend('X (Vicon)', 'Y (Vicon)', 'Z (Vicon)', 'X (Camera)', 'Y (Camera)', 'Z (Camera)', 'location','east');
xlabel('Time (s)', 'Fontsize',16);
ylabel('Position (mm)', 'Fontsize',16);
axis tight;
ax = axis;
axis([ax(1:2) ax(3)-100 ax(4)+100]);
subplot(212);
plot(v(:,1) - v(1,1), v(:,5:7), m(:,1) - v(1,1) + off(1), m(:,5:7));
legend('Roll (Vicon)', 'Pitch (Vicon)', 'Yaw (Vicon)', 'Roll (Camera)', 'Pitch (Camera)', 'Yaw (Camera)', 'location','east');
xlabel('Time (s)', 'Fontsize',16);
ylabel('Euler angles (rad)', 'Fontsize',16);
axis tight;
ax = axis;
axis([ax(1:2) ax(3)-.5 ax(4)+.5]);
print -dpdf icra17_vicon_pnp.pdf;

%% feature vector figure

dosub = false;
if dosub
    clf
end

% sample data
clear d set;
d1 = data38('abs');
d1.vel = pose_to_vel(d1.vei, d1.iws);
d1.a = 47200;
d1.b = d1.a+1500;
d{1} = d1;
d2 = data38('vinyl');
d2.vel = pose_to_vel(d2.vei, d2.iws);
d2.a = 55200;
d2.b = d2.a+1500;
d{2} = d2;

for i=1:2
    d{i}.pre = romano_features('pre', d{i}.iws, d{i}.vei, d{i}.ai, mass, 150, 0, [d{i}.a-1515 d{i}.b+1515]);
    d{i}.pre = d{i}.pre(11:20,:);
    d{i}.post = romano_features('post', d{i}.pre, 5, 'naive', 0, 1);
    d{i}.post(:, (end-3):end) = d{i}.post(:,[end-1 end end-3 end-2]); % swap V and Ft
    
    d{i}.m = mean(d{i}.post);
    d{i}.r = range(bsxfun(@minus, d{i}.post, d{i}.m));
    d{i}.m(1:5) = mean(d{i}.m(1:5));
    d{i}.r(1:5) = mean(d{i}.r(1:5))*3;
    d{i}.m([6 8]) = mean(d{i}.m([6 8]));
    d{i}.r([6 8]) = mean(d{i}.r([6 8]));
end

meanm = mean([d{1}.m; d{2}.m]);
meanr = mean([d{1}.r; d{2}.r]);
for i=1:2
    d{i}.posted = bsxfun(@rdivide, bsxfun(@minus, d{i}.post, meanm), meanr);
    minmin(i) = min(min(d{i}.posted));
    maxmax(i) = max(max(d{i}.posted));
end
minmin = min(minmin);
maxmax = max(maxmax);

f = figure;
imagesc([d{1}.posted d{2}.posted]);
ax = gca;
clim = ax.CLim;
close(f);

for i=1:2
    if dosub
        subplot(2,2, 1 + (i-1));
    else
        figure;
    end
    set(gca, 'FontSize',12);
    hold on;
    for j=0:150:1500
        line(d{i}.iws(d{i}.a+[j j],1)-d{i}.iws(1,1), [-15 30], 'color',[.6 .6 .6])
        if j < 1500
            text(d{i}.iws(d{i}.a+j+75)-d{i}.iws(1,1), -14, sprintf('%d', j/150+1), 'color',[.6 .6 .6], 'horizontalalignment','center');
        end
    end
    datums = plot(d{i}.iws(d{i}.a:d{i}.b,1)-d{i}.iws(1,1), ...
         [...
          d{i}.iws(d{i}.a:d{i}.b,4) ...
          sqrt(sum(d{i}.iws(d{i}.a:d{i}.b,2:3).^2,2)) ...
          dft321(d{i}.ai(d{i}.a:d{i}.b,2:4))*10 ...
          sqrt(sum(d{i}.vel(d{i}.a:d{i}.b,:).^2,2))/10 ...
         ]);
    axis([d{i}.iws(d{i}.a,1)-d{i}.iws(1,1) d{i}.iws(d{i}.b,1)-d{i}.iws(1,1) -15 30]);
    set(gca, 'XTick', round(d{i}.iws(d{i}.a + (225:450:1499)',1) - d{i}.iws(1,1), 2));
    legend(datums([1 2 4 3]), 'Normal force (N)', 'Tangential force (N)', 'Tip speed (cm/s)', 'Acceleration (cm/s^2)', 'location','northeast');
    xlabel('Time (s)', 'FontSize',16);
    box on;
    if ~dosub
        print('-dpdf', sprintf('icra17_fv_%d_1.pdf', i));
    end
    
    if dosub
        sub = subplot(2,2, 3 + (i-1));
    else
        figure;
        sub = axes;
    end
    
    posted = d{i}.posted;
    posted = posted';
    
    cla;
    imagesc(posted, clim);
    set(gca, 'FontSize',12);
    colormap jet;
    box off;
    sub.XTick = 1:10;
    sub.YTickLabel = [];
    xlabel('Feature vectors', 'FontSize',16);
    sub.XRuler.Axle.Visible = 'off';
    sub.YRuler.Axle.Visible = 'off';
    labels = repmat({''}, [1 5]);
    for j=1:5
        labels{j} = sprintf('%g-%g kHz', (j-1)*.3, j*.3);
    end
    things = {'Fn', 'Ft', 'V'};
    for thing=1:length(things)
        labels{end+1} = sprintf('Mean %s', things{thing});
        labels{end+1} = sprintf('Std %s', things{thing});
    end
    for j=1:length(labels)
        text(.3, j, labels{j}, 'FontSize',12, 'HorizontalAlignment','right');
    end
    if ~dosub
        print('-dpdf', sprintf('icra17_fv_%d_2.pdf', i));
    end
end

if ~dosub
    figure;
    set(colorbar, 'XTick', [0 1], 'XTickLabel', {'min' 'max'});
    colormap jet;
    print -dpdf icra17_fv_colorbar.pdf;
end

%% histograms

force14 = cell2mat(cellfun(@(f) f(:,3), features14(:,4), 'uniformoutput',false));
force38 = cell2mat(cellfun(@(f) f(:,3), features38(:,4), 'uniformoutput',false));

% manually picked limits
flim = [0 35];
slim = [0 200];

lbl = @(diam, val, units) {
    ylabel('Incidence');
    xlabel(sprintf('%s (%s)', val, units));
    title(sprintf('%s (%s end-effector)', val, diam));
};

%figure;
subplot(221);
histogram(force14(force14 > flim(1) & force14 < flim(2)), 50);
lbl('6.35 mm', 'Normal force', 'N');
subplot(223);
histogram(force38(force38 > flim(1) & force38 < flim(2)), 50);
lbl('9.525 mm', 'Normal force', 'N');
subplot(222);
histogram(speed14(speed14 > slim(1) & speed14 < slim(2)), 50);
lbl('6.35 mm', 'Tip speed', 'mm/s');
subplot(224);
histogram(speed38(speed38 > slim(1) & speed38 < slim(2)), 50);
lbl('9.525 mm', 'Tip speed', 'mm/s');

%print -dpdf icra17_histograms.pdf
