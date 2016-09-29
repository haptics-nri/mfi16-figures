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
           
int = process_mini40(100*ones(size(int,1),1), int, zeros(1,6));

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
    elseif str2num(episodes(ep).name) < 7 % change this to select end-effector
        continue;
    end
    
    flow = parse_flow([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep 'stickcam.flow']);
    material = flow.answers('surface name').text;
    fprintf('Loading %s (%s)...\n', episodes(ep).name, material);
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
bfeatures = cell(0, 5); % cols: label, vibration, speed, normal, tangential
for m = 1:length(materials)
    ep = data(materials{m});
    
    fprintf('Romano features for %s\n', materials{m});
    %%
    new_feats = romano_features('pre', ep.iws, ep.vei, ep.ai, mass, 150, [5 .5], ep.ss);
    %%
    features = [features
                num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
            
    new_feats = romano_features('pre', ep.biws, ep.bvei, ep.bai, mass, 150, [5 .5], ep.bss);
    bfeatures = [bfeatures
                num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
end

% test/train split

% 4/5 train, 1/5 test
split_idx = randsample(1:2, size(features,1), true, [4/5 1/5]);
bsplit_idx = randsample(1:2, size(bfeatures,1), true, [4/5 1/5]);

%% vicon features

train_features = features(split_idx==1, :);
test_features  = features(split_idx==2, :);

%% bluefox features

train_features = bfeatures(bsplit_idx==1, :);
test_features  = bfeatures(bsplit_idx==2, :);

%%
% crossval

cv = cvpartition(cell2mat(train_features(:,1)), 'KFold', 3);
oc_confusion = cell(1, cv.NumTestSets);
mc_confusion = cell(1, cv.NumTestSets);
oc_answers = cell(1, cv.NumTestSets);
mc_answers = cell(1, cv.NumTestSets);
%%
% hyperparameters
nbins = 20:10:60; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.4; % 25
nu = .05:0.05:0.3; % .6
gamma = 10:20:100; % 200
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
        trainrange = range(train_vectors(:,2:end));
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

%% test on test set (first run GS/CV to completion)

[~, gsi] = max(gs_acc);
gs_nbins = nbins(gs_idx(gsi,1));
gs_binmode = binmode{gs_idx(gsi,2)};
gs_alpha = alpha(gs_idx(gsi,3));
gs_nu = nu(gs_idx(gsi,4));
gs_gamma = gamma(gs_idx(gsi,5));
gs_stmode = stmode(gs_idx(gsi,6));
%%
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
trainrange = range(train_vectors(:,2:end));
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


%% feature vectors -- first set gsi to optimal and run the grid search iter
fv1 = figure;
%fv2 = figure;
f = romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode);
f = bsxfun(@minus, f, mean(f));
f = bsxfun(@rdivide, f, range(f));
allmin = min(min(f));
allmax = max(max(f));
names = {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'};
for i=1:5
    idx = cell2mat(train_features(:,1))==i;
    g = f(idx,:);
    g = [g mean(g(:,[end-5 end-3 end-1]), 2)];
    g = sortrows(g, size(g,2));
    g = g(:,1:end-1);
    g(:, (end-3):end) = g(:,[end-1 end end-3 end-2]); % swap V and Ft

    figure(fv1);
    subplot(5,1,i);
    imagesc(g, [allmin allmax]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box off; axis off;
    text(0.3, size(g,1)/2, ...
         names{i}, ...
         'FontSize', 14, ...
         'HorizontalAlignment', 'right', ...
         'Interpreter', 'tex');
     
    %figure(fv2);
    %subplot(1,6,i);
    %cor = corrcoef(g(:,1:end-6));
    %imagesc(cor);
end
figure(fv1);
a = subplot(5,1,5);
axis on;
a.XRuler.Axle.Visible = 'off';
a.YRuler.Axle.Visible = 'off';
%a.XTick = 1:10;
labels = {};
for i=1:gs_nbins
    labels{end+1} = sprintf('Bin %d', i);
end
things = {'Fn', 'Ft', 'V'};
for thing=1:length(things)
    labels{end+1} = sprintf('Mean %s', things{thing});
    if gs_stmode
        labels{end+1} = sprintf('Std %s', things{thing});
    end
end
%a.XTickLabels = labels;
%a.XTickLabelRotation = 70;
%a.TickLength = [0 0];
for i=1:length(labels)
    text(i, size(g,1)*1.2, labels{i}, 'Rotation',-90, 'FontSize',14);
end
colormap jet;
axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
set(colorbar('ticks',[]), 'edgecolor','none');
text(1.08, 1, '1', 'fontsize',14);
text(1.08, 0, '0', 'fontsize',14);
print -dpdf icra17_feature_vectors.pdf;
%figure(fv2);
%subplot(1,6,6);
%colormap jet;
%colorbar;

%% confusion matrices -- first set gsi to optimal and run the test set
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 2)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'};
ax.YTickLabel = {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'};
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
print -dpdf icra17_confusion_precision.pdf;

fprintf('\n');
fprintf('Surface & Precision & Recall & $F_1$ score\n');
for i=1:length(materials)
    fprintf('%s  &  ', material_names{i});
    others = setdiff(1:5, i);
    tp = mc_test_confusion(i,i);
    fp = sum(mc_test_confusion(others,i));
    tn = sum(sum(mc_test_confusion(others,others)));
    fn = sum(mc_test_confusion(i,others));
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    f1 = 2*prec*rec/(prec+rec);
    fprintf('%.3f  &  %.3f  &  %.3f', prec, rec, f1);
    fprintf(' \\\\ \n');
end

%% accelerometer comparison figure

[v,f,da,dg,~,a,~,dt,~,~,m] = load_stick('../../nri/data/20160906/stickcam/6/');
t = cumsum(bitand(dt, 65535))/1e6;
figure
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(211)
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(212)
plot(t, bsxfun(@minus, a(:,5:7), mean(a(:,5:7))))
subplot(211)
axis([9 17 -8 8])
subplot(212)
axis([9 17 -8 8])
xlabel('Time (s)')
ylabel('Internal accelerometer signal (m/s^2)')
subplot(211)
xlabel('Time (s)')
ylabel('External accelerometer signal (m/s^2)')
print -dpdf -r0 icra17_accel_compare.pdf

%% feature vector figure

dosub = false;
if dosub
    clf
end

% sample data
d1 = data38('abs');
d1.vel = pose_to_vel(d1.vei, d1.iws);
d1.a = 47200;
d1.b = d1.a+1500;
d{1} = d1;
d2 = data38('vinyl');
d2.vel = pose_to_vel(d2.vei, d2.iws);
d2.a = 55200;
d2.b = d2.a+1500;
d{2} = d2;

for i=1:2
    d{i}.pre = romano_features('pre', d{i}.iws, d{i}.vei, d{i}.ai, mass, 150, 0, [d{i}.a-1515 d{i}.b+1515]);
    d{i}.pre = d{i}.pre(11:20,:);
    d{i}.post = romano_features('post', d{i}.pre, 5, 'naive', 0, 1);
    d{i}.post(:, (end-3):end) = d{i}.post(:,[end-1 end end-3 end-2]); % swap V and Ft
    
    d{i}.m = mean(d{i}.post);
    d{i}.r = range(bsxfun(@minus, d{i}.post, d{i}.m));
    d{i}.m(1:5) = mean(d{i}.m(1:5));
    d{i}.r(1:5) = mean(d{i}.r(1:5))*3;
    d{i}.m([6 8]) = mean(d{i}.m([6 8]));
    d{i}.r([6 8]) = mean(d{i}.r([6 8]));
end

meanm = mean([d{1}.m; d{2}.m]);
meanr = mean([d{1}.r; d{2}.r]);
for i=1:2
    d{i}.posted = bsxfun(@rdivide, bsxfun(@minus, d{i}.post, meanm), meanr);
    minmin(i) = min(min(d{i}.posted));
    maxmax(i) = max(max(d{i}.posted));
end
minmin = min(minmin);
maxmax = max(maxmax);

f = figure;
imagesc([d{1}.posted d{2}.posted]);
ax = gca;
clim = ax.CLim;
close(f);

for i=1:2
    if dosub
        subplot(2,2, 1 + (i-1));
    else
        figure;
    end
    set(gca, 'FontSize',12);
    hold on;
    for j=0:150:1500
        line(d{i}.iws(d{i}.a+[j j],1)-d{i}.iws(1,1), [-15 30], 'color',[.6 .6 .6])
        if j < 1500
            text(d{i}.iws(d{i}.a+j+75)-d{i}.iws(1,1), -14, sprintf('%d', j/150+1), 'color',[.6 .6 .6], 'horizontalalignment','center');
        end
    end
    datums = plot(d{i}.iws(d{i}.a:d{i}.b,1)-d{i}.iws(1,1), ...
         [...
          d{i}.iws(d{i}.a:d{i}.b,4) ...
          sqrt(sum(d{i}.iws(d{i}.a:d{i}.b,2:3).^2,2)) ...
          dft321(d{i}.ai(d{i}.a:d{i}.b,2:4))*10 ...
          sqrt(sum(d{i}.vel(d{i}.a:d{i}.b,:).^2,2))/10 ...
         ]);
    axis([d{i}.iws(d{i}.a,1)-d{i}.iws(1,1) d{i}.iws(d{i}.b,1)-d{i}.iws(1,1) -15 30]);
    set(gca, 'XTick', round(d{i}.iws(d{i}.a + (225:450:1499)',1) - d{i}.iws(1,1), 2));
    legend(datums([1 2 4 3]), 'Normal force (N)', 'Tangential force (N)', 'Tip speed (cm/s)', 'Acceleration (cm/s^2)', 'location','northeast');
    xlabel('Time (s)');
    box on;
    if ~dosub
        print('-dpdf', sprintf('icra17_fv_%d_1.pdf', i));
    end
    
    if dosub
        sub = subplot(2,2, 3 + (i-1));
    else
        figure;
        sub = axes;
    end
    
    posted = d{i}.posted;
    posted = posted';
    %posted = posted - minmin;
    %posted = posted / (maxmax - minmin);
    
    cla;
    imagesc(posted, clim);
    set(gca, 'FontSize',12);
    colormap jet;
    %colorbar;
    box off;
    sub.XTick = 1:10;
    sub.YTickLabel = [];
    xlabel('Feature vectors');
    sub.XRuler.Axle.Visible = 'off';
    sub.YRuler.Axle.Visible = 'off';
    labels = repmat({''}, [1 5]);
    for j=1:5
        labels{j} = sprintf('%g-%g kHz', (j-1)*.3, j*.3);
    end
    things = {'Fn', 'Ft', 'V'};
    for thing=1:length(things)
        labels{end+1} = sprintf('Mean %s', things{thing});
        labels{end+1} = sprintf('Std %s', things{thing});
    end
    for j=1:length(labels)
        text(.3, j, labels{j}, 'FontSize',12, 'HorizontalAlignment','right');
    end
    if ~dosub
        print('-dpdf', sprintf('icra17_fv_%d_2.pdf', i));
    end
end

if ~dosub
    figure;
    set(colorbar, 'XTick', [0 1], 'XTickLabel', {'min' 'max'});
    colormap jet;
    print -dpdf icra17_fv_colorbar.pdf;
end

%% run both grid searches (cellsplit expand first)

% vicon
icra17_figures_21;
icra17_figures_23;
icra17_figures_24;
icra17_figures_29;
vgs_acc = gs_acc;

% bluefox
icra17_figures_22;
icra17_figures_23;
icra17_figures_24;
icra17_figures_29;
bgs_acc = gs_acc;

%% GS sensitivity analysis

gs_plots(vgs_acc, gs_idx(:,[1 3 4 5]), {'nbins', nbins; 'alpha', alpha; 'nu', nu; 'gamma', gamma});
suplabel('Vicon', 't');

gs_plots(bgs_acc, gs_idx(:,[1 3 4 5]), {'nbins', nbins; 'alpha', alpha; 'nu', nu; 'gamma', gamma});
suplabel('Bluefox', 't');
