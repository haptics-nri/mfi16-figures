% this file loads and verifies the data from our first two sphere
% calibration sessions

addpath(genpath('RANSAC-Toolbox'))

%% load data

% load translation and rotation from socket1 and socket4
v1 = csvload('../../nri/data/20151010/calibration.socket1/socket1.csv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
v4 = csvload('../../nri/data/20151010/calibration.socket4/socket4.csv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
         
vf = csvload('../../nri/data/20151222/socket1stick/vicon.tsv', ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});

% load Vicon model
pts = vst2points('../../nri/data/20151010/proton.vsk');

% Stu's data from CAD
stu = [  88.9     0     -18.45   % Vicon reflector positions in Vicon frame (mm)
        -71.84  -71.84  -50.2
          0    -135.26 -113.7
          0      97.16  -81.95
          0       0       0   ];
dstu = [194.47    0    -582.65]'; % end-effector pos in Vicon frame (mm)

% internal mini40
int1 = csvload('../../nri/data/20151010/calibration.socket1/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
int4 = csvload('../../nri/data/20151010/calibration.socket4/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
intf = csvload('../../nri/data/20151218/socket1stick/teensy.ft.csv', ...
               [{'Timestamp'}, ...
                arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);


%% verify Vicon model

[R, T, A, B, order, fit] = vicon_match(stu, pts);
assert(all(all(order == [1 2 3 4 5; 1 2 3 4 5])));

fprintf('Vicon reflector positions match in order with %g-deg rotation and %g-mm translation, avg. distance %g mm between matched points\n', norm(xfconv(R))*180/pi, norm(T), fit);

%% preprocess data

% narrow to actual data (magic numbers picked by eye)
x1 = v1(500:1500, 2:4);
x4 = v4(1000:6000, 2:4);
xf = vf(500:8000, 2:4);

int1 = process_mini40(int1);
int4 = process_mini40(int4);
intf = process_mini40(intf);

%% sphere fits (translation matching)

[c1, r1, e1, cs1] = sphereFit_ransac(x1, 10);
[c4, r4, e4, cs4] = sphereFit_ransac(x4, 10);
[cf, rf, ef, csf] = sphereFit_ransac(xf, 10);

%% rotation matching



%% save fits

save socket.mat v1 v4 pts stu dstu R T A B order fit x1 x4 c1 r1 e1 cs1 c4 r4 e4 cs4;

%% check goodness of fits and agreement between them

fprintf('Sphere fit to socket1: c=%s mm, r=%g mm, e=%g, %d inliers\n', mat2str(c1, 5), r1, e1, nnz(cs1));
fprintf('Sphere fit to socket4: c=%s mm, r=%g mm, e=%g, %d inliers\n', mat2str(c4, 5), r4, e4, nnz(cs4));
fprintf('Distance between centers = %g mm\n', norm(c1 - c4));
fprintf('Difference in radius = %g mm\n', abs(r1 - r4));

proj1 = zeros(nnz(cs1), 3);
j = 0;
for i=find(cs1)
    j = j + 1;
    proj1(j,:) = R\(xfconv(v1(500+i-1,5:7))\(c1 - v1(500+i-1,2:4))' - T);
end
avg1 = mean(proj1)';
proj4 = zeros(nnz(cs4), 3);
j = 0;
for i=find(cs4)
    j = j + 1;
    proj4(j,:) = R\(xfconv(v4(1000+i-1,5:7))\(c4 - v4(1000+i-1,2:4))' - T);
end
avg4 = mean(proj4)';
fprintf('Averages:\n');
fprintf('\tsocket1=%s (avg std %.2f) (%.2f mm, %.2f deg from dstu)\n', mat2str(avg1, 5), mean(std(proj1)), norm(avg1 - dstu), acos(dot(avg1, dstu)/norm(avg1)/norm(dstu))*180/pi);
fprintf('\tsocket4=%s (avg std %.2f) (%.2f mm, %.2f deg from dstu)\n', mat2str(avg4, 5), mean(std(proj4)), norm(avg4 - dstu), acos(dot(avg4, dstu)/norm(avg4)/norm(dstu))*180/pi);
fprintf('\tdstu   =%s\n', mat2str(dstu, 5));
fprintf('\tinter  = %.2f mm, %.2f deg\n', norm(avg1 - avg4), acos(dot(avg1, avg4)/norm(avg1)/norm(avg4))*180/pi);

clf;
hold on;
sphereplot(c1, r1);
sphereplot(c4, r4);
hold off;
