% this file loads and verifies the data from the sphere-recalibration

addpath(genpath('RANSAC-Toolbox'))

%% load data

% vicon data
v1 = csvload('../../nri/data/20160223/socket1stick/vicon.tsv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
v2 = csvload('../../nri/data/20160223/socket2stick/vicon.tsv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
v3 = csvload('../../nri/data/20160223/socket3stick/vicon.tsv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
         
% vicon model
pts = vst2points('../../nri/data/20151010/proton.vsk');

% Stu's data from CAD
stu = [  88.9     0     -18.45   % Vicon reflector positions in Vicon frame (mm)
        -71.84  -71.84  -50.2
          0    -135.26 -113.7
          0      97.16  -81.95
          0       0       0   ];
dstu = [194.47    0    -582.65]'; % end-effector pos in Vicon frame (mm)

% internal mini40
int1 = csvload('../../nri/data/20160223/socket1stick/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
int2 = csvload('../../nri/data/20160223/socket2stick/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
int3 = csvload('../../nri/data/20160223/socket3stick/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);

%% verify Vicon model

[R, T, A, B, order, fit] = vicon_match(stu, pts);
assert(all(all(order == [1 2 3 4 5; 1 2 3 4 5])));

fprintf('Vicon reflector positions match in order with %g-deg rotation and %g-mm translation, avg. distance %g mm between matched points\n', norm(xfconv(R))*180/pi, norm(T), fit);

%% preprocess data

% use the second and third (ball popped up in first)
x = [v2; v3];

int = [process_mini40(int2); process_mini40(int3)];

%% sphere fits (translation matching)

[c, r, e, cs] = sphereFit_ransac(x(:,2:4), 50); % FIXME allowing a huge amount of noise to get a reasonable number of inliers

%% check goodness of fit and plot

fprintf('Sphere fit: c=%s mm, r=%g mm, e=%g, %d inliers\n', mat2str(c, 5), r, e, nnz(cs));

clf;
hold on;
sphereplot(c, r, {x(cs,2:4) x(~cs,2:4)});
hold off;
grid on; axis equal vis3d;
