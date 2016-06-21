%% load data

%datadir = '../../../data';
%epdate = '20150804';
%epname = 'spherecalib2';
datadir = '../../nri/data';
epdate = '20151103';
epname = 'free1opto';
%%
warning('off', 'MATLAB:table:ModifiedVarnames');

fprintf('Loading episode "%s" from %s...', epname, datestr(datenum(epdate, 'yyyymmdd'), 'mm/dd/yyyy'));
fprintf('(vicon)...');
%matlab_parser_is_dumb = strsplit(epname, '.');
%vicon = csvload([datadir filesep epdate filesep epname filesep matlab_parser_is_dumb{end} '.csv'], ...
vicon = csvload([datadir filesep epdate filesep epname filesep 'vicon.tsv'], ...
                {'Timestamp', ...
                 'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
                 'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
                {'Delimiter', '\t'});
fprintf('(acc)...');
acc  = csvload([datadir filesep epdate filesep epname filesep 'teensy.acc.csv' ], ...
               {'Timestamp', 'FIFOPosition', 'AccX', 'AccY', 'AccZ'});
fprintf('(gyro)...');
gyro = csvload([datadir filesep epdate filesep epname filesep 'teensy.gyro.csv'], ...
               {'Timestamp', 'FIFOPosition', 'GyroX', 'GyroY', 'GyroZ'});
fprintf('(mag)...');
mag  = csvload([datadir filesep epdate filesep epname filesep 'teensy.mag.csv' ], ...
               {'Timestamp', 'MagX', 'MagY', 'MagZ'});
%%
fprintf('(structure)...');
stmeta = csvload([datadir filesep epdate filesep epname filesep 'structure_times.csv'], ...
                 {'UnixTimestamp', 'FrameNumber'});
             % stmeta(:,1) = timestamps, stmeta(:,2) = frame numbers
structure = zeros(size(stmeta,1), 480, 640); % structure = row x col x frame number
for i=1:size(stmeta,1)
    structure(i,:,:) = imread([datadir filesep epdate filesep epname filesep sprintf('structure%d.png', stmeta(i,2))]);
end
fprintf('done.\n');

%%

save([epdate '-' epname '.mat'], '-v7.3', 'acc', 'gyro', 'mag', 'stmeta', 'structure');

%% interpolate IMU data

acc  = unfifo(acc);
gyro = unfifo(gyro);

%% ICP

% from http://forums.structure.io/t/getting-colored-point-cloud-data-from-aligned-frames/4094/2 with probable typo (in _fy) fixed
fx = 305.73/320*640;
fy = 305.62/240*480;
cx = 159.69/320*640;
cy = 119.86/240*480;

% allocate memory
ps = cell(size(structure,1),1);
fs = cell(size(structure,1),1);
vs = cell(size(structure,1),1);
rs = zeros(size(structure,1)-1, 3,3);
ts = zeros(size(structure,1)-1, 3);
pairs = cell(size(structure,1)-1,1);

ps{1} = detectSURFFeatures(squeeze(structure(1,:,:)));
[fs{1}, vs{1}] = extractFeatures(squeeze(structure(1,:,:)), ps{1});

for i=2:size(structure,1)
    disp(i);
    
    ps{i} = detectMSERFeatures(squeeze(structure(i,:,:)));
    [fs{i}, vs{i}] = extractFeatures(squeeze(structure(i,:,:)), ps{i});
    
    pairs{i-1} = matchFeatures(fs{i-1}, fs{i});
    m1 = round(vs{i-1}(pairs{i-1}(:,1), :).Location);
    m2 = round(vs{i}(pairs{i-1}(:,2), :).Location);
    
    xyz1 = zeros(size(pairs{i-1},1), 3);
    xyz2 = zeros(size(pairs{i-1},1), 3);
    for j=1:size(pairs{i-1},1)
        xyz1(j,1) = structure(i-1,m1(j,2),m1(j,1))/10000 * (m1(j,2) - cx)/fx;
        xyz1(j,2) = structure(i-1,m1(j,2),m1(j,1))/10000 * (cy - m1(j,1))/fy;
        xyz1(j,3) = structure(i-1,m1(j,2),m1(j,1))/10000;
        
        xyz2(j,1) = structure(i,m2(j,2),m2(j,1))/10000 * (m2(j,2) - cx)/fx;
        xyz2(j,2) = structure(i,m2(j,2),m2(j,1))/10000 * (cy - m2(j,1))/fy;
        xyz2(j,3) = structure(i,m2(j,2),m2(j,1))/10000;
    end
    
    [rs(i-1,:,:), ts(i-1,:)] = rigid_transform_3D(xyz1, xyz2);
end

%%
xs = zeros(size(ts));
for i=2:size(xs,1)
    xs(i,:) = squeeze(rs(i-1,:,:))*xs(i-1,:)' + ts(i-1,:)';
end

%% EKF from A Tale of Two Helicopters

% You'll notice this section is empty.
% That's because there is precious little detail given in the paper, and
% they undermine the meager explanation at the end by saying "...actually we
% didn't implement it this way." What a waste of a publication.

%% EKF from MEAM 620

% "model B: quadrotor with good acceleration sensor"
% state: position, euler angles,    velocity,   gyro bias,     acc bias
%        x, y, z,  phi, theta, psi, vx, vy, vz, bgx, bgy, bgz, bax, bay, baz
ng = [0 0 0]';
na = [0 0 0]';
nbg = [0 0 0]';
nba = [0 0 0]';
R = @(q) [cos(q(3))*cos(q(2)) - sin(q(1))*sin(q(3))*sin(q(2))   -cos(q(1))*sin(q(3))    cos(q(3))*sin(q(2)) + cos(q(2))*sin(q(1))*sin(q(3))
          cos(q(2))*sin(q(3)) + cos(q(3))*sin(q(1))*sin(q(2))    cos(q(1))*sin(q(3))    sin(q(3))*sin(q(2)) - cos(q(3))*cos(q(2))*sin(q(1))
          -cos(q(1))*sin(q(2))                                   sin(q(1))              cos(q(1))*cos(q(2))                                ];
G = @(q) [cos(q(2))     0     -cos(q(1))*sin(q(2))
          0             1      sin(q(1))
          sin(q(2))     0      cos(q(1))*cos(q(2))];
g = [0 0 -9.8]';
ff = @(a, w, p, q, v, bg, ba) [v; G(q)\(w - bg - ng); g + R(q)*(a - ba - na); nbg; nba];
f = @(a, w) @(x) ff(a, w, x(1:3), x(4:6), x(7:9), x(10:12), x(13:15));
hh = @(p, q, v, bg, ba) [p; q; v];
h = @(x) hh(x(1:3), x(4:6), x(7:9), x(10:12), x(13:15));
Q = eye(15);
R = eye(9);
x = zeros(15, size(acc,1));
P = zeros(15, 15, size(acc,1));
P(:,:,1) = eye(15);

for i=2:size(gyro,1)
    if mod(i, 1000) == 0
        fprintf('%d\n', i);
    end
    gyro_ind = i;
    [~, acc_ind] = min(abs(acc(:,1) - gyro(i,1)));
    % FIXME provide input!!
    [x(:,i), P(:,:,i)] = ekf(f(acc(acc_ind,2:4)', gyro(gyro_ind,2:4)'), x(:,i-1), P(:,:,i-1), h, zeros(9,1), Q, R);
    P(:,:,i) = nearestSPD(P(:,:,i));
end

%% UKF from ESE 650

% state: position, linear velocity, orientation, angular velocity
%        x/y/z     x/y/z            quaternion   x/y/z

x = zeros(13, size(gyro,1));
P = zeros(13, 13, size(gyro,1));
x(:,1) = [ [0 0 0]'; [0 0 0]'; [1 0 0 0]'; [0 0 0]' ];
P(:,:,1) = eye(13);

Q = 0.1*eye(13);
R = 0.01*eye(6);

vec2quat = @(v) [cos(norm(v)/2)
                 v/((norm(v)==0) + norm(v))*sin(norm(v)/2)];
quatmul = @(q, r) [r(1)*q(1) - r(2)*q(2) - r(3)*q(3) - r(4)*q(4)
                   r(1)*q(2) + r(2)*q(1) - r(3)*q(4) + r(4)*q(3)
                   r(1)*q(3) + r(2)*q(4) + r(3)*q(1) - r(4)*q(2)
                   r(1)*q(4) - r(2)*q(3) + r(3)*q(2) + r(4)*q(1)];

ff = @(dt, p, v, q, w) [p + v*dt
                        v
                        quatmul(vec2quat(w*dt), q)
                        w];
f = @(dt) @(x) ff(dt, x(1:3), x(4:6), x(7:10), x(11:13));

hh = @(dt, p, v, q, w) [v/dt
                        w];
h = @(dt) @(x) hh(dt, x(1:3), x(4:6), x(7:10), x(11:13));

accfac = 9.8/mean(sqrt(sum(acc(:,2:4).^2,2)));
gyrofac = 0.002;

for i=2:size(gyro,1)
    if mod(i,1000) == 0
        disp(i);
    end
    
    gi = i;
    [~, ai] = min(abs(acc(:,1) - gyro(i,1)));
    
    dt = gyro(gi,1) - gyro(gi-1,1);
    z = [accfac*acc(ai,2:4)'; gyrofac*gyro(gi,2:4)'];
    
    [x(:,i), P(:,:,i)] = ukf(f(dt), x(:,i-1), P(:,:,i-1), h(dt), z, Q, R);
    %x(7:10,i) = x(7:10,i)/norm(x(7:10,i));
end
