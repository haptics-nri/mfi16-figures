%% generates figures for the ICRA 2017 paper

% calibrates using the datasets taken on 2/26/16, 7/26/16, and 8/11/16
% loads the dataset taken on 8/15/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

%DATADIR = '/Volumes/shared/Projects/Proton Pack/Data';
DATADIR = '../../nri/data';

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
[~,f] = load_stick([DATADIR filesep date filesep 'weigh/1/']);
[~,r] = sphereFit_ransac(f(:,2:4)); % 98.9% inliers

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

% dataset parameters
date = '20160906';
flowtype = 'stickcam';

% end-effector mass comes from calibration above

data = containers.Map;
episodes = dir([DATADIR filesep date filesep flowtype]);
for ep = 1:length(episodes)
    if episodes(ep).name(1) == '.'
        continue;
    end
    
    flow = parse_flow([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep 'stickcam.flow']);
    material = flow.answers('surface name').text;
    [v, int, da,dg,mi, acc, ma, dt, o,b, motrak] = ...
        load_stick([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep]);
    data(material) = struct('v',v, 'int',int, 'acc',acc, ...
                            'imu', struct('acc',da, 'gyro',dg, 'mag',ma), ...
                            'sensor', struct('opto',o, 'bio',b), ...
                            'motrak',motrak, 'dt',dt);
end

materials = data.keys;

if false
    for m = 1:length(materials)
        % replace NUC timestamps with Teensy deltas
        ep = data(materials{m});
        ep.nuc_t = ep.int(:,1);
        ep.int(:,1) = ep.int(1,1) + cumsum(bitand(ep.dt, 65535))/1e6;
        data(materials{m}) = ep;
    end
end

%% process with vicon

for m = 1:length(materials)
    ep = data(materials{m});
    
    [~,~,~,~,~,~, vei, ai, ~,~,~, iws] = process_stick(ep.v, ep.int, ep.acc, mass, [0;0;0], H_vic2bod, H_m402bod, H_bal2imu, -ep.off);
    ep.vei = vei;
    ep.ai = ai;
    ep.iws = iws;
    
    data(materials{m}) = ep;
end

%% process with bluefox

for m = 1:length(materials)
    ep = data(materials{m});
    
    [~,~,~,~,~,~, vei, ai, ~,~,~, iws] = process_stick(ep.motrak, ep.int, ep.acc, mass, [0;0;0], H_vic2bod, H_m402bod, H_bal2imu);
    ep.bvei = vei;
    ep.bai = ai;
    ep.biws = iws;
    
    data(materials{m}) = ep;
end

%% SVM stuff (following Romano & KJK 2014 + Strese & Schuwerk & Steinbach 2015)

% DO NOT RUN THIS AGAIN -- TRAINING ON TEST SET WILL RESULT

% extract features
features = cell(0, 5); % cols: label, vibration, speed, normal, tangential
for m = 1:length(materials)
    ep = data(materials{m});
    
    fprintf('Romano features for %s\n', materials{m});
    %%
    new_feats = romano_features('pre', ep.iws, ep.vei, ep.ai, mass, 0.05, [20 3]);
                                                                 % FIXME reexamine these thresholds
    %%
    features = [features
                num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
end

% test/train split

% 4/5 train, 1/5 test
split_idx = randsample(1:2, size(features,1), true, [4/5 1/5]);

train_features = features(split_idx==1, :);
test_features  = features(split_idx==2, :);

%%
% crossval

cv = cvpartition(cell2mat(train_features(:,1)), 'KFold', 3);
oc_confusion = cell(1, cv.NumTestSets);
mc_confusion = cell(1, cv.NumTestSets);
oc_answers = cell(1, cv.NumTestSets);
mc_answers = cell(1, cv.NumTestSets);
%%
% hyperparameters
nbins = 1:2:10; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.5; % 25
nu = .05:0.05:0.3; % .6
gamma = 1:2:20; % 200
stmode = true; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for i=2:size(gs_idx,1)
    gs_idx(i,:) = gs_idx(i-1,:);
    for j=size(gs_idx,2):-1:1
        if gs_idx(i,j) == gs_limits(j)
            gs_idx(i,j) = 1;
        else
            gs_idx(i,j) = gs_idx(i,j) + 1;
            break;
        end
    end
end
%%
gs_acc = zeros(size(gs_idx,1),1);
clear romano_features; % clear persistent vars
elapsed = tic;
for gsi=1:size(gs_idx,1)
    gs_nbins = nbins(gs_idx(gsi,1));
    gs_binmode = binmode{gs_idx(gsi,2)};
    gs_alpha = alpha(gs_idx(gsi,3));
    gs_nu = nu(gs_idx(gsi,4));
    gs_gamma = gamma(gs_idx(gsi,5));
    gs_stmode = stmode(gs_idx(gsi,6));
    %%
    iter = tic;
    fprintf('Grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g, stmode=%d\n', gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma, gs_stmode);
    cv_acc = zeros(cv.NumTestSets,1);
    for cvi = 1:cv.NumTestSets
        train_vectors = [cell2mat(train_features(cv.training(cvi),1)) ...
                         romano_features('post', train_features(cv.training(cvi),2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];
        val_vectors   = [cell2mat(train_features(cv.test(cvi),    1)) ...
                         romano_features('post', train_features(cv.test(cvi)    ,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];

        trainmean = mean(train_vectors(:,2:end));
        train_vectors(:,2:end) = bsxfun(@minus, ...
                                        train_vectors(:,2:end), ...
                                        trainmean);
        val_vectors  (:,2:end) = bsxfun(@minus, ...
                                        val_vectors  (:,2:end), ...
                                        trainmean);
        trainrange = max(train_vectors(:,2:end)) - min(train_vectors(:,2:end));
        train_vectors(:,2:end) = bsxfun(@rdivide, ...
                                        train_vectors(:,2:end), ...
                                        trainrange);
        val_vectors  (:,2:end) = bsxfun(@rdivide, ...
                                        val_vectors  (:,2:end), ...
                                        trainrange);

        models = cell(5, ... % vibration, speed+normal, normal+tangential, speed+normal+tangential, all
                      length(materials)+1); % OC models for each materials, one MC model
        common_args = ' -q ';
        %oc_train_args = ['-s 2 -t 2 -n 0.0303' common_args];
        mc_train_args = [sprintf('-m 1000 -s 1 -t 2 -n %g -g %g', gs_nu, gs_gamma) common_args];
        test_args = common_args;

        %for mi=1:length(materials)
        %    fprintf('Material: %s\n', materials{mi});
        %    %%
        %    % normalize features
        %    train_labels =  ones(nnz(train_vectors(:,1)==mi), 1);
        %    val_labels   =  ones(nnz(val_vectors  (:,1)==mi), 1);
        %    unval_labels = -ones(nnz(val_vectors  (:,1)~=mi), 1);
        %    train_feats = train_vectors(train_vectors(:,1)==mi, 2:end);
        %    val_feats   = val_vectors  (val_vectors  (:,1)==mi, 2:end);
        %    unval_feats = val_vectors  (val_vectors  (:,1)~=mi, 2:end);
        %
        %    % train SVM
        %    gamma = 0.0303;%evangelista(train_feats);
        %    models{mi}                    = svmtrain(  train_labels, train_feats, sprintf('%s -g %g', train_args, gamma));
        %    [in_pred, in_acc, in_prob]    = svmpredict(val_labels  , val_feats  , models{mi}, test_args);
        %    [out_pred, out_acc, out_prob] = svmpredict(unval_labels, unval_feats, models{mi}, test_args);
        %    fprintf('\tin-class accuracy: %g%%\n' , 100*nnz(in_pred  == 1)/length(in_pred));
        %    fprintf('\tout-class accuracy: %g%%\n', 100*nnz(out_pred ~= 1)/length(out_pred));
        %end

        models{end} = svmtrain(train_vectors(:,1), train_vectors(:,2:end), mc_train_args);

        if ~isempty(models{end})
            % evaluate by comparing all OCSVMs and the MCSVM
            %oc_confusion{cvi} = zeros(length(materials));
            mc_confusion{cvi} = zeros(length(materials));
            %prob = zeros(size(val_vectors,1),length(materials));
            %for mi=1:length(materials)
            %    prob(:,mi) = rabaoui_dissim(models{mi}, val_vectors(:,2:end));
            %end
            %[~, oc_answers{cvi}] = min(prob, [], 2);
            mc_answers{cvi} = svmpredict(zeros(size(val_vectors,1),1), val_vectors(:,2:end), models{end}, '-q');

            for i=1:length(materials)
                for j=1:length(materials)
            %        oc_confusion{cvi}(i,j) = nnz(oc_answers{cvi}(val_vectors(:,1)==i) == j);
                    mc_confusion{cvi}(i,j) = nnz(mc_answers{cvi}(val_vectors(:,1)==i) == j);
                end
            end

            cv_acc(cvi) = sum(diag(mc_confusion{cvi}))/sum(sum(mc_confusion{cvi}));
            fprintf('\tFold %d: MC %g%%\n', cvi, 100*cv_acc(cvi));
        else
            fprintf('\tFold %d: failed to train\n', cvi);
            cv_acc(cvi) = 0;
        end
    end
    gs_acc(gsi) = mean(cv_acc);
    fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi), toc(iter), toc(elapsed));
end
