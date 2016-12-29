%% setup

addpath(genpath('libsvm'));
addpath('../steinbach/htk-mfcc');
DATADIR = '../steinbach/LMT_TextureDB';
NAMES = readtable(fullfile(DATADIR, 'names.csv'));
materials = table2cell(NAMES(:,8));
files = dir(fullfile(DATADIR, 'Training', 'Accel', [materials{1} '_query*.txt']));
train_ids = cellfun(@(s) subsref(regexp(s, '_query(.*)\.txt', 'tokens'), substruct('{}', {1})), {files.name});
files = dir(fullfile(DATADIR, 'Testing', 'AccelScansDFT321', [materials{1} '_*.txt']));
test_ids = cellfun(@(s) subsref(regexp(s, '_(.*)\.txt', 'tokens'), substruct('{}', {1})), {files.name});

%% load data

train_episodes = struct('material', [], ...
                        'id',       [], ...
                        'acc',      [], ...
                        'fric',     [], ...
                        'snd',      []);
test_episodes = struct('material', [], ...
                       'id',       [], ...
                       'acc',      [], ...
                       'fric',     [], ...
                       'snd',      []);

m = 0;
while m < length(materials)
    m = m + 1;
    for i=1:length(train_ids)
        t = tic;
        fprintf('Loading training data for %d/%d %s (%d/%d)... ', m, length(materials), materials{m}, i, length(train_ids));
        try
            [acc, fric, snd] = load_lmt(DATADIR, 'train', materials{m}, train_ids{i});
        catch
            fprintf('skip!\n');
            materials(m) = [];
            m = m - 1;
            break;
        end
        train_episodes((m-1)*length(train_ids) + i) = struct('material', materials{m}, ...
                                                             'id',       train_ids{i}, ...
                                                             'acc',      acc, ...
                                                             'fric',     fric, ...
                                                             'snd',      snd);
        fprintf('%g s\n', toc(t));
    end
end
%save -v7.3 st_train.mat train_episodes;
%clear train_episodes;

if false
    for m=1:length(materials)
        for i=1:length(test_ids)
            fprintf('Loading testing data for %d/%d %s (%d/%d)...\n', m, length(materials), materials{m}, i, length(test_ids));
            [acc, fric, snd] = load_lmt(DATADIR, 'test', materials{m}, test_ids{i});
            test_episodes((m-1)*length(test_ids) + i) = struct('material', materials{m}, ...
                                                               'id',       test_ids{i}, ...
                                                               'acc',      acc, ...
                                                               'fric',     fric, ...
                                                               'snd',      snd);
        end
    end
    %save -v7.3 st_test.mat test_episodes;
    %clear test_episodes;
end

%% train

feats = {'MF', 'H', 'SC', 'TR', 'SR', 'WV', 'SP', 'F', 'RG', 'Fr', 'FM', 'SIH', 'SILH', 'SISR', 'SISH', 'SISS'};

if ~exist('train_episodes', 'var')
    load st_train;
end

for i=1:length(train_episodes)
    fprintf('Calculating features for training episode %d/%d... ', i, length(train_episodes));
    t = tic;
    train_episodes(i).features = steinbach_features(feats, train_episodes(i));
    fprintf('%g s\n', toc(t));
end

%% test

%% divide training set
train_features = [];
val_features = [];
train_labels = [];
val_labels = [];
for i = 1:length(train_episodes)/10
    d = randsample(10, 7);
    ti = (i-1)*10 + d;
    vi = (i-1)*10 + setdiff(1:10, d);
    train_features = [train_features; cell2mat({train_episodes(ti).features}')];
    val_features = [val_features; cell2mat({train_episodes(vi).features}')];
    train_labels = [train_labels cellfun(@(m) find(strcmp(m, materials), 1), {train_episodes(ti).material})];
    val_labels = [val_labels cellfun(@(m) find(strcmp(m, materials), 1), {train_episodes(vi).material})];
end

%% I dunno try a fucking SVM I guess
train_features(isnan(train_features)) = 0;
val_features(isnan(val_features)) = 0;
trainmean = mean(train_features);
train_vectors = bsxfun(@minus, ...
                       train_features, ...
                       trainmean);
val_vectors   = bsxfun(@minus, ...
                       val_features  , ...
                       trainmean);
trainrange = range(train_vectors);
train_vectors = bsxfun(@rdivide, ...
                       train_vectors, ...
                       trainrange);
val_vectors   = bsxfun(@rdivide, ...
                       val_vectors  , ...
                       trainrange);

model = svmtrain(train_labels', train_vectors, '-m 1000 -s 1 -t 2 -n 0.1 -g 0.5 -q');
answers = svmpredict(zeros(size(val_labels')), val_vectors, model, '-q');
confusion = zeros(length(materials));
for i=1:length(materials)
    for j=1:length(materials)
        confusion(i,j) = nnz(answers(val_labels == i) == j);
    end
end
fprintf('Test set accuracy: %g\n', sum(diag(confusion))/sum(sum(confusion)));

