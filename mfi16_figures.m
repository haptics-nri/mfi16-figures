%% generates tbl:regression for the MFI 2016 paper

% calibrates using the datasets taken on 2/26/16 and 2/23/16
% loads the datasets taken on 2/19/16 and 3/10/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

DATADIR = '../../nri/data';

%% sphere calibration (see go_sphere_again.m)

% vicon data
v2 = csvload([DATADIR '/20160223/socket2stick/vicon.tsv'], ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
v3 = csvload([DATADIR '/20160223/socket3stick/vicon.tsv'], ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});

% use the second and third (ball popped up in first)
x = [v2; v3];

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

date = '20160226'; %'20151218';
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
freecalib.int_s = int_s;
freecalib.grav = grav;
freecalib.ideal_grav = ideal_grav;
freecalib.c = c;
freecalib.r = r;
freecalib.cs = inliers;
freecalib.R = R;
freecalib.t = t;

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
materials = {'black', 'white', 'blue', 'brown', 'red'};
material_names = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
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

v =    cell(length(materials), length(reps), length(tools));
int  = cell(length(materials), length(reps), length(tools));
acc  = cell(length(materials), length(reps), length(tools));
gyr  = cell(length(materials), length(reps), length(tools));
mic  = cell(length(materials), length(reps), length(tools));
aacc = cell(length(materials), length(reps), length(tools));
%%
for mi = 1:length(materials)
    for ri = 1:length(reps)
        for ti = 1:length(tools)
            %%
            fprintf('Loading data for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            dataset = [materials{mi} reps{ri} tools{ti}];
            prefix = [DATADIR filesep date{ri} filesep dataset filesep];
            
            [v{mi,ri,ti}, int{mi,ri,ti}, acc{mi,ri,ti}, gyr{mi,ri,ti}, mic{mi,ri,ti}, aacc{mi,ri,ti}] = load_stick(prefix);
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
                = process_stick(v{mi,ri,ti}, int{mi,ri,ti}, aacc{mi,ri,ti}, mass, com, H_vic2bod, H_m402bod, H_bal2imu, offset(ri));
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

% DO NOT RUN THIS AGAIN -- TRAINING ON TEST SET WILL RESULT
% instead load featuresplit20160329.mat

% extract features
features = cell(0, 5); % cols: label, vibration, speed, normal, tangential
for mi = 1:length(materials)
    for ri = 1:length(reps)/2
        for ti = 1:length(tools)
            fprintf('Romano features for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            %%
            new_feats = romano_features('pre', intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, accint{mi,ri,ti}, mass, 0.05, [20 3]);
                                                                         % FIXME reexamine these thresholds
            %%
            features = [features
                        num2cell(repmat(mi, size(new_feats,1), 1)) new_feats];
        end
    end
end

% test/train split

% 4/5 train, 1/5 test
split_idx = randsample(1:2, size(features,1), true, [4/5 1/5]);

train_features = features(split_idx==1, :);
test_features  = features(split_idx==2, :);

%%
% crossval

cv = cvpartition(cell2mat(train_features(:,1)), 'KFold', 5);
oc_confusion = cell(1, cv.NumTestSets);
mc_confusion = cell(1, cv.NumTestSets);
oc_answers = cell(1, cv.NumTestSets);
mc_answers = cell(1, cv.NumTestSets);
%%
% hyperparameters
nbins = 1:2:20; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.5; % 25
nu = .05:0.05:0.7; % .6
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

%%
clf;
subplot(1,2,1);
bar3(mc_confusion);
title(sprintf('MC accuracy = %g%%', 100*sum(diag(mc_confusion))/sum(sum(mc_confusion))));
subplot(1,2,2);
bar3(oc_confusion);
title(sprintf('OC accuracy = %g%%', 100*sum(diag(oc_confusion))/sum(sum(oc_confusion))));

%% test on test set (first run GS/CV to completion)

[~, gsi] = max(gs_acc);
gs_nbins = nbins(gs_idx(gsi,1));
gs_binmode = binmode{gs_idx(gsi,2)};
gs_alpha = alpha(gs_idx(gsi,3));
gs_nu = nu(gs_idx(gsi,4));
gs_gamma = gamma(gs_idx(gsi,5));
gs_stmode = stmode(gs_idx(gsi,6));

train_vectors = [cell2mat(train_features(:,1)) ...
                          romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];
test_vectors  = [cell2mat(test_features(:,1)) ...
                          romano_features('post', test_features (:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];
trainmean = mean(train_vectors(:,2:end));
train_vectors(:,2:end) = bsxfun(@minus, ...
                                train_vectors(:,2:end), ...
                                trainmean);
test_vectors (:,2:end) = bsxfun(@minus, ...
                                test_vectors (:,2:end), ...
                                trainmean);
trainrange = max(train_vectors(:,2:end)) - min(train_vectors(:,2:end));
train_vectors(:,2:end) = bsxfun(@rdivide, ...
                                train_vectors(:,2:end), ...
                                trainrange);
test_vectors (:,2:end) = bsxfun(@rdivide, ...
                                test_vectors (:,2:end), ...
                                trainrange);
mc_train_args = [sprintf('-m 1000 -s 1 -t 2 -n %g -g %g', gs_nu, gs_gamma) common_args];
final_model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), mc_train_args);
mc_test_answers = svmpredict(zeros(size(test_vectors,1),1), test_vectors(:,2:end), final_model, test_args);
mc_test_confusion = zeros(length(materials));
for i=1:length(materials)
    for j=1:length(materials)
        mc_test_confusion(i,j) = nnz(mc_test_answers(test_vectors(:,1)==i) == j);
    end
end
fprintf('Test set accuracy: %g\n', sum(diag(mc_test_confusion))/sum(sum(mc_test_confusion)));

%% save

save mfi16_tbl_regression.mat

%% figures

clear set; % there is a var called "set" but I need to use the set() function

% sphere calibration
figure;
sphereplot(spherecalib.c, spherecalib.r, {spherecalib.x(:,2:4)});
xlabel('X position (mm)')
ylabel('Y position (mm)')
zlabel('Z position (mm)')
print -dpdf mfi16_sphere_calib.pdf

% free space calibration
figure;
plot(freecalib.int_s(:,1)-freecalib.int_s(1,1), freecalib.int_s(:,2:4))
xlabel('Time (s)')
ylabel('Measured force (N)')
legend('X component', 'Y component', 'Z component')
print -dpdf mfi16_freespace_data.pdf
figure;
transformed_grav = freecalib.ideal_grav;
for i=1:size(freecalib.ideal_grav,1)
    transformed_grav(i,2:4) = freecalib.R * freecalib.ideal_grav(i,2:4)' + freecalib.t;
end
hold on;
h1 = plot3(freecalib.grav(:,2), freecalib.grav(:,3), freecalib.grav(:,4), '.');
h2 = plot3(transformed_grav(:,2), transformed_grav(:,3), transformed_grav(:,4), '.');
grid on
axis equal vis3d
view(11.5, 42)
xlabel X
ylabel Y
zlabel Z
legend('Measured gravity in body frame', 'World-frame gravity transformed to body frame', 'location','southeast')
print -dpdf mfi16_freespace_grav.pdf

% typical dataset
figure;
subplot(2,1,1);
plot(vbodyint{3,2,1}(:,1)-vbodyint{3,2,1}(1,1), vbodyint{3,2,1}(:,2:4));
yl = ylabel('Position (mm)');
yl.Position(1) = yl.Position(1) - .3;
legend('X', 'Y', 'Z', 'location','east');
set(gca, 'fontsize',14);
subplot(2,1,2);
plot(intworldsub{3,2,1}(:,1)-intworldsub{3,2,1}(1,1), intworldsub{3,2,1}(:,2:4));
xlabel('Time (s)')
yl = ylabel('Force (N)');
yl.Position(1) = yl.Position(1) - .3;
legend('X', 'Y', 'Z', 'location','northeast');
set(gca, 'fontsize',14);
print -dpdf mfi16_typical_data.pdf

%% feature vectors -- first set gsi to optimal and run the grid search iter
fv1 = figure;
fv2 = figure;
f = romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode);
m = mean(f);
for i=1:5
    figure(fv1);
    subplot(5,1,i);
    g = f(cell2mat(train_features(:,1))==i, :);
    g = bsxfun(@minus, g, m);
    g = bsxfun(@rdivide, g, max(g) - min(g));
    g = [g mean(g(:,(end-6):end), 2)];
    g = sortrows(g, size(g,2));
    g = g(:,1:end-1);
    imagesc(g);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box off; axis off;
    text(0.3, size(g,1)/2, ...
         material_names{i}, ...
         'FontSize', 14, ...
         'HorizontalAlignment', 'right', ...
         'Interpreter', 'tex');
     
    figure(fv2);
    subplot(1,6,i);
    cor = corrcoef(g(:,1:end-6));
    imagesc(cor);
end
figure(fv1);
colormap jet;
axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
set(colorbar('ticks',[]), 'edgecolor','none');
print -dpdf mfi16_feature_vectors.pdf;
figure(fv2);
subplot(1,6,6);
colormap jet;
colorbar;

%% confusion matrices -- first set gsi to optimal and run the test set
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 1)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.YTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.FontSize = 14;
xlabel('Detected material');
ylabel('Actual material');
ax.XLabel.Position = ax.XLabel.Position + [0 0.1 0];
for i=1:length(materials)
    for j=1:length(materials)
        if i == j
            c = 'white';
        else
            c = 'black';
        end
        text(j, i, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(:,j))), ...
             'FontSize',14, 'Color',c, 'horizontalalignment','center');
    end
end
print -dpdf mfi16_confusion_precision.pdf;
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 2)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.YTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.FontSize = 14;
xlabel('Detected material');
ylabel('Actual material');
ax.XLabel.Position = ax.XLabel.Position + [0 0.1 0];
for i=1:length(materials)
    for j=1:length(materials)
        if i == j
            c = 'white';
        else
            c = 'black';
        end
        text(j, i, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(i,:))), ...
             'FontSize',14, 'Color',c, 'HorizontalAlignment','center');
    end
end
print -dpdf mfi16_confusion_recall.pdf;

fprintf('\n');
fprintf('Surface & Accuracy & Precision & Recall & $F_1$ score\n');
for i=1:length(materials)
    fprintf('%s  &  ', material_names{i});
    others = setdiff(1:5, i);
    tp = mc_test_confusion(i,i);
    fp = sum(mc_test_confusion(others,i));
    tn = sum(sum(mc_test_confusion(others,others)));
    fn = sum(mc_test_confusion(i,others));
    acc = (tp+tn)/(tp+tn+fp+fn);
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    f1 = 2*prec*rec/(prec+rec);
    fprintf('%.3f  &  %.3f  &  %.3f  &  %.3f', acc, prec, rec, f1);
    fprintf(' \\\\ \n');
end
