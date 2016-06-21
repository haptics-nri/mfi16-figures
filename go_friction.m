%% setup

addpath(genpath('RANSAC-Toolbox'))
addpath(genpath('geom3d'))

%% load data

dir = '../../nri/data';
date = '20160113';
material = 'table';
tool = 'opto';
rep = '1';

% end-effector properties
mass = 0.1094; % kg
com = [-.056; .014; .0218]; % m

dataset = [material rep tool];
prefix = [dir filesep date filesep dataset filesep];

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

opto = csvload([prefix 'optoforce.csv'], ...
               {'Timestamp', 'X', 'Y', 'Z'});



%% process data

int = process_mini40(int);
acc = unfifo(acc);


%% visualize vicon and force

%%% NB we are using the IMU frame as the body frame!

H_vic2acc = [ 1       0       0         5.65
              0      -1       0         6.57
              0       0      -1      -559.01
              0       0       0         1   ];
H_imu2m40 = [ 0       0       1       108.99
             -1       0       0         0.53
              0      -1       0        -2.98
              0       0       0         1   ];
H_imu2bal = [ 1       0       0       188.83
              0       1       0         0.57
              0       0       1        23.63
              0       0       0         1   ];
H_fit =     [ 1       0       0       246.96
              0       1       0       -23.33
              0       0       1       -13.9
              0       0       0         1   ];
H_imu2opf = [ 0       0.9063  0.4226  191.96
             -1       0       0         0.57
              0      -0.4226  0.9063   30.36
              0       0       0         1   ];
H_imu2bio = [ 0.9397  0      -0.342   191.08
              0      -1        0        0.57
             -0.342   0       -0.9397  24.04
              0       0        0        1   ];
   
% transform Vicon into body frame and end-effector frame
vbody = v;
vend = v;
for i=1:size(vbody,1)
    tf = xfconv(v(i,2:7)) * H_vic2acc;
    vbody(i,2:4) = tf(1:3,4);
    vbody(i,5:7) = xfconv(tf(1:3,1:3));
    tf = tf * H_fit;                       %%% FIXME just use the tooling ball vector
    vend(i,2:4) = tf(1:3,4);
    vend(i,5:7) = xfconv(tf(1:3,1:3));
end

%%%%% FIXME Vicon-Mini40 time offset found by hand! what's wrong with the clock synchronization???
offset = 25.0047;

% upsample and offset Vicon/accel to match Mini40
vint     = [int(:,1) interp1(v(:,1),     v(:,2:4),     int(:,1)-offset) slerp(v(:,1),     v(:,5:7),     int(:,1)-offset)];
vbodyint = [int(:,1) interp1(vbody(:,1), vbody(:,2:4), int(:,1)-offset) slerp(vbody(:,1), vbody(:,5:7), int(:,1)-offset)];
vendint  = [int(:,1) interp1(vend(:,1),  vend(:,2:4),  int(:,1)-offset) slerp(vend(:,1),  vend(:,5:7),  int(:,1)-offset)];
accint   = [int(:,1) interp1(acc(:,1),   acc(:,2:4),   int(:,1)-offset)];
%%
% transform Mini40 into body frame and world frame
intbody = int;
intworld = int;
for i=1:size(intbody,1)
    intbody(i,2:4) = H_imu2m40(1:3,1:3)*int(i,2:4)';
    %intbody(i,5:7) = H_imu2m40(1:3,1:3)*int(i,5:7)';
    tf = xfconv(vbodyint(i,5:7));
    intworld(i,2:4) = tf * intbody(i,2:4)';
    %intworld(i,5:7) = tf * intbody(i,5:7)';
end
%%
% subtract weight of end-effector
intworldsub = intbody;
for i=1:size(intbody,1)
    fg = mass * [0; 0; -9.81];
    intworldsub(i,2:4) = intworld(i,2:4) - fg';
    %intworldsub(i,5:7) = intworld(i,5:7) - cross(com, fg)';
end
%%
% this plots the total trajectory with an arrow representing the current
% position + orientation of the rig, and a slider to move back and forth
clf;
plot3(v(:,2), v(:,3), v(:,4), vbody(:,2), vbody(:,3), vbody(:,4), vint(:,2), vint(:,3), vint(:,4), vbodyint(:,2), vbodyint(:,3), vbodyint(:,4));
hold on;
quiver3(vbodyint(1:100:end,2), vbodyint(1:100:end,3), vbodyint(1:100:end,4), intworldsub(1:100:end,2), intworldsub(1:100:end,3), intworldsub(1:100:end,4), 5);
R = xfconv(vint(1,5:7));
Rbody = xfconv(vbodyint(1,5:7));
hat = quiver3([vint(1,2) vbodyint(1,2)], ...
              [vint(1,3) vbodyint(1,3)], ...
              [vint(1,4) vbodyint(1,4)], ...
              [R(1,3) -Rbody(1,3)], ...
              [R(2,3) -Rbody(2,3)], ...
              [R(3,3) -Rbody(3,3)], ...
              0.2, 'linewidth',3);
vec = plot3([vint(1,2) vbodyint(1,2)], ...
            [vint(1,3) vbodyint(1,3)], ...
            [vint(1,4) vbodyint(1,4)], ...
            'linewidth', 2);
legend('', '', 'Vicon frame', 'Forces', 'Mini40 frame', 'Up vector', 'Proton pack', 'location','best');
hold off; grid on; axis equal vis3d; view(30, 30);
q = @(i,R,Rbody) {set(hat, 'XData',[vint(i,2) vbodyint(i,2)], ...
                           'YData',[vint(i,3) vbodyint(i,3)], ...
                           'ZData',[vint(i,4) vbodyint(i,4)], ...
                           'UData',[R(1,3) -Rbody(1,3)], ...
                           'VData',[R(2,3) -Rbody(2,3)], ...
                           'WData',[R(3,3) -Rbody(3,3)])
                  set(vec, 'XData',[vint(i,2) vbodyint(i,2)], ...
                           'YData',[vint(i,3) vbodyint(i,3)], ...
                           'ZData',[vint(i,4) vbodyint(i,4)])};
p = @(s,e) {q(round(s.Value), xfconv(vint(round(s.Value),5:7)), xfconv(vbodyint(round(s.Value),5:7)))
            title(sprintf('t = %g s', vint(round(s.Value),1)-vint(1,1)))};
uicontrol('Style','slider', ...
          'Min',1, 'Max',size(intbody,1), 'Value',1, ...
          'Position',[100 10 400 20], ...
          'CreateFcn',p, 'Callback',p)
uicontrol('Style','text', ...
          'Position',[75 10 20 20], ...
          'String','0');
uicontrol('Style','text', ...
          'Position',[505 10 50 20], ...
          'String',sprintf('%g',vint(end,1)-vint(1,1)));
      
%% pretty pictures

      
%% filtering and stuff

a = .01;
b = [1 .01-1];

vendfilt = filtfilt(a, b, vendint(:,2:4));
speed = bsxfun(@rdivide, diff(vendfilt), diff(vendint(:,1)));
speed = [speed; 0 0 0]; % make it the same length
speedfilt = filtfilt(a, b, speed);

accfilt = filtfilt(a, b, acc(:,2:4));

forcemask = sqrt(sum(intworldsub(:,2:4).^2,2)) > 2;
speedmask = sqrt(sum(speedfilt.^2,2)) > 20;

maskplane = fitPlane(vendint(forcemask,2:4));

fn = zeros(1,size(speed,1));
ft = zeros(1,size(speed,1));
for i=1:size(intworldsub,1)
    fn(i) = abs(dot(intworldsub(i,2:4)-.1953, planeNormal(maskplane)));
    ft(i) = abs(dot(intworldsub(i,2:4)-.1953, speedfilt(i,:)/norm(speedfilt(i,:))));
            %%% FIXME project into plane and then onto velocity
end
fn = filtfilt(a, b, fn);
ft = filtfilt(a, b, ft);
      
%% okay let's just measure some frequencies



%%

%[x] plot force rotated into world frame (not IMU frame)
%    - switch H_imu2m40: intbody = H * int
%    - rotate to world: intworld = vbody \ intbody

%[x] fit plane to end-effector position (not IMU position)
%    - use a mask to fit to only the time with nonzero forces

%[x] measure speed in end-effector frame (not IMU frame)
%[x] low pass filter
%[x] choose threshold
%[ ] static: plot Fp vs Ft, should go thru origin
%[ ] kinetic: plot Fp vs Ft, should be flat
%[x] project onto velocity direction

%[x] maybe interp1 the magnitude + direction of rotations separately? (use SLERP)




%[x] show magnitude of fit error in gravity compensation -- is it really a perfect sphere?
%[x] figure out why plot3(vendint(forcemask, 2:4)) doesn't look right
%[?] fix offset in intworldsub



%GFTT
%- friction coefficient
%- spectral centroid of tap transient



% fix the stupid frames
% make a fancy animation like heather's