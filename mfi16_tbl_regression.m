%% generates tbl:regression for the MFI 2016 paper

% calibrates using the datasets taken on 2/26/16 and 2/23/16
% loads the datasets taken on 2/19/16 and 3/10/16 and does machine learning

%% sphere calibration (see go_sphere_again.m)

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

% vicon data
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

% use the second and third (ball popped up in first)
x = [v2; v3];

[c, ~, ~, cs] = sphereFit_ransac(x(:,2:4), 50); % FIXME allowing a huge amount of noise to get a reasonable number of inliers

d = nan([size(x,1) 3]);
for i = find(cs)
    cc = xfconv(x(i,2:7)) \ [c 1]';
    d(i,:) = cc(1:3);
end
d = nanmean(d);

% the product of sphere calibration is d

%% free-space calibration (see go_bias.m)

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
           
int = process_mini40(int, zeros(1,6));

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
H1 = [R [0 y z]';      0 0 0 1];
H2 = [eye(3) [x 0 0]'; 0 0 0 1];
S = solve(H1 * H2 * [0 0 0 1]' == [d 1]');

H_vic2bod(1:3,1:3) = R';
H_vic2bod(2,4) = S.y;
H_vic2bod(3,4) = S.z;
H_bal2imu(1,4) = S.x;

%% setup for learning

% dataset parameters
dir = '../../nri/data';
materials = {'black', 'white', 'blue', 'brown', 'red'};
tools = {'stick'};
reps = {'1', '2', '3', '4', '5', '1', '2', '3', '4', '5'};
date = [repmat({'20160310'}, 1, 5) repmat({'20160219'}, 1, 5)];

%                 reps
% 2016-03-10 2016-02-19
video_offsets = [20  0 13 10  9 24 44 23 32 25
                 17 34 25 27 27 42 22 36 44 36 % materials
                 23 25 18 21 24 37 29 50 58 34
                 26 20 22 21 31 26 29 50 53 48
                 26 19 27 22 29 61 46 37 32 46];

% end-effector properties
mass = 0.1503; % kg
com = [-.0029; -.00301; .0348]; % m
%%

v =   cell(length(materials), length(reps), length(tools));
int = cell(length(materials), length(reps), length(tools));
acc = cell(length(materials), length(reps), length(tools));
%%
for mi = 1:length(materials)
    for ri = 1:length(reps)
        for ti = 1:length(tools)
            %%
            fprintf('Loading data for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            dataset = [materials{mi} reps{ri} tools{ti}];
            prefix = [dir filesep date{ri} filesep dataset filesep];
            
            [v{mi,ri,ti}, int{mi,ri,ti}, acc{mi,ri,ti}] = load_stick(prefix);
        end
    end
end

%% preprocessing

offset = [repmat(-2.5887, 1, 5) repmat(5.9088, 1, 5)]; % 2016-03-10 2016-02-19
%%
vv =          cell(length(materials), length(reps), length(tools));
ii =          cell(length(materials), length(reps), length(tools));
vbody =       cell(length(materials), length(reps), length(tools));
vend =        cell(length(materials), length(reps), length(tools));
vint =        cell(length(materials), length(reps), length(tools));
vbodyint =    cell(length(materials), length(reps), length(tools));
vendint =     cell(length(materials), length(reps), length(tools));
intbody =     cell(length(materials), length(reps), length(tools));
intworld =    cell(length(materials), length(reps), length(tools));
intworldsub = cell(length(materials), length(reps), length(tools));
accint =      cell(length(materials), length(reps), length(tools));
accworld =    cell(length(materials), length(reps), length(tools));
vel =         cell(length(materials), length(reps), length(tools));
%%
for mi = 1:length(materials)
    for ri = 1:length(reps)
        for ti = 1:length(tools)
            %%
            fprintf('Processing data for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            [vv{mi,ri,ti}, ii{mi,ri,ti}, ...
             vbody{mi,ri,ti}, vend{mi,ri,ti}, ...
             vint{mi,ri,ti}, vbodyint{mi,ri,ti}, vendint{mi,ri,ti}, ...
             accint{mi,ri,ti}, accworld{mi,ri,ti}, ...
             intbody{mi,ri,ti}, intworld{mi,ri,ti}, intworldsub{mi,ri,ti}] ...
                = process_stick(v{mi,ri,ti}, int{mi,ri,ti}, acc{mi,ri,ti}, mass, com, H_vic2bod, H_m402bod, H_bal2imu, offset(ri));
        end
    end
end

%% features

mu_k   = cell(length(materials), length(reps), length(tools));
spring = cell(length(materials), length(reps), length(tools));
power  = cell(length(materials), length(reps), length(tools));
%%
for mi = 1:length(materials)
    for ri = 1:length(reps)
        for ti = 1:length(tools)
            %%
            fprintf('Extracting features for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            fprintf('\tfriction coefficient\n');
            [mu, err] = extract_friction(intworldsub{mi,ri,ti}, vendint{mi,ri,ti});
            mu_k{mi,ri,ti} = real([mu err]);
            
            fprintf('\tspringiness\n');
            [k, z0, err] = extract_springiness(intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, mass);
            spring{mi,ri,ti} = [k z0 err];
            
            fprintf('\tvibration power\n');
            power{mi,ri,ti} = extract_vibration(intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, accint{mi,ri,ti}, mass);
        end
    end
end

%% SVM stuff (following Romano & KJK 2014)

% extract features
features = zeros(0, 34); % first col is labels
for mi = 1:length(materials)
    for ri = 1:length(reps)
        for ti = 1:length(tools)
            fprintf('Romano features for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            %%
            new_feats = romano_features(intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, accint{mi,ri,ti}, mass, 0.1, [20 5]);
            %%
            features = [features
                        repmat(mi, size(new_feats,1), 1) new_feats];
        end
    end
end

% test/train split

% 3/5 train, 1/5 validation, 1/5 test
split_idx = randsample(1:3, size(features,1), true, [3/5 1/5 1/5]);

train_features = features(split_idx==1, :);
val_features   = features(split_idx==2, :);
test_features  = features(split_idx==3, :);

trainmean = mean(train_features(:,2:end));
trainvar  = var (train_features(:,2:end));
train_features(:,2:end) = bsxfun(@rdivide, ...
                                 bsxfun(@minus, ...
                                        train_features(:,2:end), ...
                                        trainmean), ...
                                 trainvar);
val_features  (:,2:end) = bsxfun(@rdivide, ...
                                 bsxfun(@minus, ...
                                        val_features  (:,2:end), ...
                                        trainmean), ...
                                 trainvar);
test_features (:,2:end) = bsxfun(@rdivide, ...
                                 bsxfun(@minus, ...
                                        test_features (:,2:end), ...
                                        trainmean), ...
                                 trainvar);
%%
% ML
models   = cell(1, length(materials)+1);
common_args = ' -q ';
train_args = ['-s 2 -t 2 -n 0.8' common_args];
test_args = common_args;

for mi=1:length(materials)
    fprintf('Material: %s\n', materials{mi});
    %%
    % normalize features
    train_labels =  ones(nnz(train_features(:,1)==mi), 1);
    val_labels   =  ones(nnz(val_features  (:,1)==mi), 1);
    unval_labels = -ones(nnz(val_features  (:,1)~=mi), 1);
    train_feats = train_features(train_features(:,1)==mi, 2:end);
    val_feats   = val_features  (val_features  (:,1)==mi, 2:end);
    unval_feats = val_features  (val_features  (:,1)~=mi, 2:end);
    
    % train SVM
    gamma = 0.0303;%evangelista(train_feats);
    models{mi}                    = svmtrain(  train_labels, train_feats, sprintf('%s -g %g', train_args, gamma));
    [in_pred, in_acc, in_prob]    = svmpredict(val_labels  , val_feats  , models{mi}, test_args);
    [out_pred, out_acc, out_prob] = svmpredict(unval_labels, unval_feats, models{mi}, test_args);
    fprintf('\tin-class accuracy: %g%%\n' , 100*nnz(in_pred  == 1)/length(in_pred));
    fprintf('\tout-class accuracy: %g%%\n', 100*nnz(out_pred ~= 1)/length(out_pred));
end

models{end} = svmtrain(train_features(:,1), train_features(:,2:end), '-s 1 -t 2 -n 0.5 -g 10 -q');

% evaluate by comparing all OCSVMs
oc_confusion = zeros(length(materials));
mc_confusion = zeros(length(materials));
en_confusion = zeros(length(materials));
for i=1:size(val_features,1)
    prob = zeros(1,length(materials));
    for mi=1:length(materials)
        prob(mi) = rabaoui_dissim(models{mi}, val_features(i,2:end));
    end
    [~, oc_answer] = min(prob);
    oc_confusion(val_features(i,1), oc_answer) = oc_confusion(val_features(i,1), oc_answer) + 1;
    
    mc_answer = svmpredict(0, val_features(i,2:end), models{end}, '-q');
    mc_confusion(val_features(i,1), mc_answer) = mc_confusion(val_features(i,1), mc_answer) + 1;
end

clf;
subplot(1,2,1);
bar3(mc_confusion);
title(sprintf('MC accuracy = %g%%', 100*sum(diag(mc_confusion))/sum(sum(mc_confusion))));
subplot(1,2,2);
bar3(oc_confusion);
title(sprintf('OC accuracy = %g%%', 100*sum(diag(oc_confusion))/sum(sum(oc_confusion))));

%% save

save mfi16_tbl_regression.mat

%% problems

% complex numbers / spikes
% 1-3
% 1-5
% 2-1
% 2-5
% 3-1
% 3-3
% 4-1
% 4-3
% 4-5
% 5-1
% 5-4

% max iter
% 1-5
% 4-3
% 4-7
