%% load data

v = csvload('../../nri/data/20151028/acrylic1opto/vicon.tsv', ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

int = csvload('../../nri/data/20151028/acrylic1opto/teensy.ft.csv', ...
              [{'Timestamp'}, ...
               arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
ext = csvread('../../nri/data/20151028/acrylic1opto/stb.csv', 1,0);
ext = ext(:,1:7);

%% process internal data

% ported from STB code
%   notes from tracing through the STB code:
%     - add sensor data to array (RunSTB.py:128)
%         - STB::read_data (hapticsstb.py:158)
%             - serial_data (hapticsstb_rt.pyx:62)
%                 - serial_m40 (hapticsstb_rt.pyx:29)
%                     - merge the bytes
%                     - handle negatives
%                     - scale down
%                     - subtract bias
%                     - transpose
%                     - multiply with transform
intM40tf =  [ 0.00679   0.01658  -0.04923   6.20566   0.15882  -6.19201 
              0.11638  -7.31729  -0.04322   3.54949  -0.08024   3.57115 
             10.35231   0.32653  10.61091   0.29668  10.33382   0.25761 
              0.00022  -0.04140   0.14917   0.02435  -0.15234   0.01567 
             -0.16837  -0.00464   0.08561  -0.03311   0.08763   0.03721 
              0.00128  -0.08962   0.00085  -0.08785   0.00204  -0.08790 ];
intbias =   [ 0         0         0         0         0         0       ]; % FIXME this is bogus

extM40tf =  [  0.165175269, 	6.193716635,	-0.05972626,	0.020033203, 	-0.136667224, 	-6.42215241	
 0.002429674, 	-3.63579423,	0.466390998, 	7.308900211, 	-0.18369186, 	-3.65179797	
 -10.5385017,	0.802731009,	-10.1357248,	0.359714766,	-10.0934065,	0.442593679	
 0.144765089,	-0.032574325,	0.004132077,	0.038285567, 	-0.145061852,	-0.010347366
 -0.089833077,	-0.024635731,	0.165602185,	-0.009131771,	-0.080132747,	0.039589968	
 0.001846317,	0.085776855,	0.005262967,	0.088317691, 	0.001450272,	0.087714269	];
extbias =   [ 0.19460667  0.10940667  0.05095667 -0.1919      0.58522    -0.13414333];

% merge bytes (preserve timestamp column), handle negatives, and scale
int = [int(:,1) bitshift(int(:,[2 4 6 8 10 12]), 8) + int(:,[3 5 7 9 11 13])];
for i=1:size(int,1)
    for j=2:size(int,2)
        if int(i,j) >= 2048
            int(i,j) = int(i,j) - 4096;
        end
    end
end
int(:,2:end) = int(:,2:end) * 0.002;
int = permute(int, [1 7 6 5 4 3 2]); % why is this necessary?

% bias and transform
int(:,2:end) = (intM40tf * bsxfun(@minus, int(:,2:end), intbias)')';

% fix conversion bug in external
ext(:,2:end) = bsxfun(@plus, (extM40tf \ ext(:,2:end)')', extbias);
for i=1:size(ext,1)
    for j=2:size(ext,2)
        if ext(i,j) >= 2048*.002
            ext(i,j) = ext(i,j) - 4096*.002;
        end
    end
end
ext(:,2:end) = (extM40tf * bsxfun(@minus, ext(:,2:end), extbias)')';

%% metric definition

metric = @(i,e) sum(sqrt(sum((i - e).^2, 2)))/size(i,1);

%% align in time

% hyperparameters
chop = 8; % seconds to chop off the end (to avoid clipped region)
skew = -0.0765; % clock skew between computers
%%
skews = -.1:.0005:.1;
back_metric = zeros(size(skews));
for i=1:length(skews)
skew = skews(i);
fprintf('%d/%d: %f\n', i, length(skews), skew);
%%
% find boundaries
% use the external Tx as a zero marker (FIXME?)
lim = mean(ext(1:10000, 5)) * 20;
a = find(abs(ext(:,5)) > lim, 1, 'first');
b = find(abs(ext(:,5)) > lim, 1, 'last');
[~, b] = min(abs(ext(:,1) - (ext(b,1)-chop)));
% chop both at the same times
[~, aa] = min(abs(int(:,1) - (ext(a,1)+skew)));
[~, bb] = min(abs(int(:,1) - (ext(b,1)+skew)));

% shrink to fit
ext_s = ext(a:b, :);
int_s = int(aa:bb, :);

% upsample external to match internal
ext_re = int_s;
ext_re(:,2:end) = interp1(ext_s(:,1), ext_s(:,2:end), int_s(:,1), 'linear', 'extrap');
ext_s = ext_re;

% remove Vicon rotation from external
for tm=1:size(ext_s,1)
    [~, tv] = min(abs(v(:,1) - ext_s(tm,1)));
    Rvw = xfconv(v(tv,5:7));
    ext_s(tm,2:4) = Rvw\ext_s(tm,2:4)';
    % FIXME deal with torque
end

%% backward fit: use forces to extract rotation

[R, t] = rigid_transform_3D(-int_s(:,2:4), ext_s(:,2:4));

back_metric(i) = metric(bsxfun(@minus, t, R*int_s(:,2:4)')', ext_s(:,2:4));

% -I = R^{-1} (E - t)
% I = -R^{-1} (E - t)
% I = R^{-1} (t - E)
% I = R^{-1} t - R^{-1} E
% I - R^{-1} t = -R^{-1} E
% R^{-1} t - I = R^{-1} E
% t - R I = E

end

%% forward fit: use rotation from Vicon to align Mini40 data

% Rvw = from Vicon frame to world frame (varies, measured by Vicon)
% Rvi = from Vicon frame to IMU frame (constant, learned by earlier calibration)
% Rio = from IMU frame to OptoForce frame (constant, given from CAD)

Hio =  [ 0  0.9063  0.4226  191.96
        -1  0.0000  0.0000    0.57
         0 -0.4226  0.9063   30.36
         0  0.0000  0.0000    1.00 ];
Rio = Hio(1:3,1:3);

Him =   [ 0  0 1 108.99
         -1  0 0    0.53
          0 -1 0   -2.98
          0  0 0     1   ];
Rim = Him(1:3,1:3);

Rvi =  [  -0.0749    0.9963   -0.0414
           0.9926    0.0784    0.0924
           0.0953   -0.0341   -0.9949 ];
%%       
% assumption: external Mini40 is aligned with world frame
% full rotation from external Mini40 to internal Mini40 is Rvw*Rvi*Rim

extrot = ext_s;
for i=1:size(extrot,1)
    tm = i;
    [~, tv] = min(abs(v(:,1) - extrot(tm,1)));
    Rvw = xfconv(v(tv,5:7));
    extrot(tm,2:4) = Rvw*Rvi*Rim*extrot(tm,2:4)';
    extrot(tm,5:7) = Rvw*Rvi*Rim*extrot(tm,5:7)';
end
