%% generates tbl:ftcal for the RSS 2016 paper

% loads the datasets taken on 12/23/15 and uses them to calculate:
%  - bias (in force and torque) of the Mini40 sensor
%  - mass of each of the four end-effectors (stick, opto, bio, center)

%% setup

addpath(genpath('RANSAC-Toolbox'))

datadir = '../../nri/data/20151223';
endeffs = {'stick', 'opto', 'bio', 'center'};

%% force calibration

mass = zeros(size(endeffs));
fbias = zeros(0, 3);

for i=1:length(endeffs)
    episodes = dir([datadir filesep 'free*' endeffs{i}]);
    for j=1:length(episodes)
        fprintf('%s\n', episodes(j).name);
        [mass_tmp, bias_tmp] = weigh([datadir filesep episodes(j).name]);
        mass(i) = mass(i) + mass_tmp;
        fbias = [fbias; bias_tmp];
    end
    mass(i) = mass(i)/length(episodes);
end
fbias_std = std(fbias);
fbias = mean(fbias);

%% torque calibration

com = zeros(length(endeffs), 3);
tbias = zeros(0, 3);

for i=1:length(endeffs)
    episodes = dir([datadir filesep 'free*' endeffs{i}]);
    for j=1:length(episodes)
        fprintf('%s\n', episodes(j).name);
        
        % TODO make this a function
        int = csvload([datadir filesep episodes(j).name filesep 'teensy.ft.csv'], ...
                      [{'Timestamp'}, ...
                       arrayfun(@(x) ['FT' num2str(x)], 0:11, 'UniformOutput',false)]);
               
        a = round(size(int,1)*1/10);
        b = round(size(int,1)*9/10);
        int = int(a:b,:);
        int = process_mini40(int, zeros(1,6));
        
        N = size(int,1);
        A = zeros(3*N, 6);
        b = zeros(3*N, 1);
        F = bsxfun(@minus, int(:,2:4), fbias); % subtract out previously measured force bias
        for k=1:N
            A( (k-1)*3+1 : k*3 , :) = [eye(3) [ 0       F(k,3) -F(k,2)
                                               -F(k,3)  0       F(k,1)
                                                F(k,2) -F(k,1)  0      ]];
            b( (k-1)*3+1 : k*3 ) = int(k,5:7)';
        end
        x = robustfit(A, b, '', '', 'off'); % robust least squares with no constant term
        
        com(i,:) = com(i,:) + x(4:6)';
        tbias = [tbias; x(1:3)'];
    end
    com(i,:) = com(i,:)/length(episodes);
end
tbias_std = std(tbias);
tbias = mean(tbias);

%% print results

fprintf('REPORT\n');
for i=1:length(endeffs)
    fprintf('%s: mass = %s g, CoM = %s mm\n', endeffs{i}, mat2str(mass(i)*1000, 4), mat2str(com(i,:)*1000, 3));
end
fprintf('Force bias: %s +/- %s N\n', mat2str(fbias, 3), mat2str(fbias_std, 3));
fprintf('Torque bias: %s +/- %s N\n', mat2str(tbias, 3), mat2str(tbias_std, 3));
fprintf('REMEMBER TO SAVE RESULTS for process.m\n');
fprintf('\t>> stamp = datetime;\n');
fprintf('\t>> save mini40_calib.mat stamp fbias tbias fbias_std tbias_std mass com endeffs;\n');
