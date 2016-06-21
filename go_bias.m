%% load data

addpath(genpath('RANSAC-Toolbox'))

dir = '../../nri/data';
date = '20160226'; %'20151218';
offset = 22.848; % manual (from taps)
material = 'free';
tool = 'stick';
rep = '1';

set = [material rep tool];
prefix = [dir filesep date filesep set filesep];

v = csvload([prefix 'vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

int = csvload([prefix 'teensy.ft.csv'], ...
              [{'Timestamp'}, ...
               arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
           
acc = csvload([prefix 'teensy.acc.csv'], ...
              {'Timestamp', 'FIFOPosition', 'AccX', 'AccY', 'AccZ'});
gyro = csvload('../../nri/data/20151104/inlab2opto/teensy.gyro.csv', ...
               {'Timestamp', 'FIFOPosition', 'GyroX', 'GyroY', 'GyroZ'});
if acc; acc = unfifo(acc); end
if gyro; gyro = unfifo(gyro); end

opto = csvload([prefix 'optoforce.csv'], ...
               {'Timestamp', 'X', 'Y', 'Z'});
           
%% process internal
int = process_mini40(int, zeros(1,6));

a = round(size(int,1)*1/5);
b = round(size(int,1)*4/5);
int_s = int(a:b,:);
if v
    [~,aa] = min(abs(v(:,1)-int(a,1)));
    [~,bb] = min(abs(v(:,1)-int(b,1)));
    v_s = v(aa:bb,:);
end
if acc
    [~,aa] = min(abs(acc(:,1)-int(a,1)));
    [~,bb] = min(abs(acc(:,1)-int(b,1)));
    acc_s = acc(aa:bb,:);
end
if gyro
    [~,aa] = min(abs(gyro(:,1)-int(a,1)));
    [~,bb] = min(abs(gyro(:,1)-int(b,1)));
    gyro_s = gyro(aa:bb,:);
end

%% find internal bias

[c,r,e,inliers] = sphereFit_ransac(int_s(:,2:4)); % 94% inliers

%% plot sphere

clf;

subplot(2,1,1);
sphereplot(c, r, {int_s(inliers,2:4), int_s(~inliers,2:4)});

subplot(2,1,2);
plot(int_s(inliers,1)-int_s(1,1), sqrt(sum(bsxfun(@minus, int_s(inliers,2:4), c).^2,2)) - r, '.', ...
     int_s(~inliers,1)-int_s(1,1), sqrt(sum(bsxfun(@minus, int_s(~inliers,2:4), c).^2,2)) - r, '.', ...
     int_s(:,1)-int_s(1,1), acos(int_s(:,3)-mean(int_s(:,3))));
%% compare with gravity

H_imu2body = [ 0       0       1       108.99
               1       0       0         0.53
               0       1       0        -2.98
               0       0       0         1   ];

[~, start] = min(abs(int_s(:,1) - v(1,1)));

grav = zeros(size(int_s,1), 4);
grav(:,1) = int_s(:,1) - offset;
ideal_grav = grav;
for tm=1:size(grav,1)
    [delta, tv] = min(abs(v(:,1) - grav(tm,1)));
    assert(delta < 0.02);
    Rvw = xfconv(v(tv,5:7));
    grav(tm,2:4) = H_imu2body(1:3,1:3) * int_s(tm,2:4)';
    ideal_grav(tm,2:4) = Rvw\[0 0 -r]';
end

[R, t, err] = rigid_ransac(ideal_grav(inliers,2:4), grav(inliers,2:4), 10); % 100% inliers??

%% TORQUE!

% 2 3    3 2
% 3 1  - 1 3
% 1 2    2 1

% f = bias + mg
% tau = bias + r x f
% bias = tau - r x f
% biasX = [tauX - (rY fZ - rZ fY)     = tauX + rZ fY
% biasY    tauY - (rZ fX - rX fZ)     = tauY - rZ fX
% biasZ    tauZ - (rX fY - rY fX)]    = tauZ
%                                     ^ if load CoM is on Z axis
%    tau is known
%    f is known
%    r is unknown (know magnitude but not direction)
%    bias is unknown
   
% |tau| = |r| |f| sin(theta)


% Goal: find bias of torque sensor
% Problem: unknown bias + unknown vector to CoM of load => too many unknowns
% Potential solution: ascertain vector to CoM of load:
%     1. "find" direction of vector by moving CoM to sensor Z axis
%     2. find magnitude of vector by solving |tau| = |r| |f| sin(theta) using force data
%               |r| = |tau|/(|f| sin(theta))
% New problem: have bias, want to find CoM of other loads
% Potential solution: { tau = bias + r x f } is three equations, three unknowns

%% find torque bias using least squares

N = size(int_s,1);
A = zeros(3*N, 6);
b = zeros(3*N, 1);
F = bsxfun(@minus, int_s(:,2:4), c); % subtract out previously measured force bias
for i=1:N
    A( (i-1)*3+1 : i*3 , :) = [eye(3) [ 0       F(i,3) -F(i,2)
                                       -F(i,3)  0       F(i,1)
                                        F(i,2) -F(i,1)  0      ]];
    b( (i-1)*3+1 : i*3 ) = int_s(i,5:7)';
end
x = robustfit(A, b, '', '', 'off'); % robust least squares with no constant term

% torque bias is x(1:3)
% vector to load CoM is x(4:6)

%% get magnitude of vector to CoM, assuming CoM is on Z axis

% |r| = |tau|/(|f| sin(theta))

f = bsxfun(@minus, int_s(:,2:4), c);
tau = bsxfun(@minus, int_s(:,5:7), x(1:3)');

dir_r = [0 0 1]; % CoM is on sensor Z axis

% cos(theta) = F `dot` r
theta = acos(f*dir_r' ./ sqrt(sum(f.^2,2)));

mag_r = sqrt(sum(tau.^2,2)) ./ (sqrt(sum(f.^2,2)) .* sin(theta));
