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

% all the features that I implemented!
feats = {'MF', 'H', 'SC', 'TR', 'SR', 'WV', 'SP', 'F', 'RG', 'Fr', 'FM', 'SIH', 'SILH', 'SISR', 'SISH', 'SISS'};

if ~exist('train_episodes', 'var')
    load st_train;
end

for i=1:length(train_episodes)
    fprintf('Calculating features for training episode %d/%d... ', i, length(train_episodes));
    t = tic;
    train_episodes(i).features = steinbach_features(feats, [],[],[], train_episodes(i));
    fprintf('%g s\n', toc(t));
end
dbclear if infnan;

%% SVM cross-validation

% construct CV splits by leaving out one recording in each split
train_features = cell(10,1);
train_labels = cell(10,1);
val_features = cell(10,1);
val_labels = cell(10,1);
% i: CV split (1..10)
% j: material (1..69)
% k: recording (1..10)
for i=1:10
    train_features{i} = [];
    val_features{i} = [];
    for j=1:length(train_episodes)/10
        for k=1:10
            idx = (j-1)*10 + k;
            if k == i % in CV split #i, leave out the i'th recording of each material
                val_features{i} = [val_features{i}; train_episodes(idx).features];
                val_labels{i} = [val_labels{i}; find(strcmp(train_episodes(idx).material, materials), 1)];
            else
                train_features{i} = [train_features{i}; train_episodes(idx).features];
                train_labels{i} = [train_labels{i}; find(strcmp(train_episodes(idx).material, materials), 1)];
            end
        end
    end
end

cvacc = zeros(1,10);
confusion = cell(1,10);
for s=1:10
    fprintf('cv %d/10\n', s);
    trainmean = mean(train_features{s});
    train_vectors = bsxfun(@minus, ...
                           train_features{s}, ...
                           trainmean);
    val_vectors   = bsxfun(@minus, ...
                           val_features{s}  , ...
                           trainmean);
    trainrange = range(train_vectors);
    train_vectors = bsxfun(@rdivide, ...
                           train_vectors, ...
                           trainrange);
    val_vectors   = bsxfun(@rdivide, ...
                           val_vectors  , ...
                           trainrange);

    model = svmtrain(train_labels{s}, train_vectors, '-m 1000 -s 1 -t 2 -n 0.1 -g 0.5 -q');
    answers = svmpredict(zeros(size(val_labels{s})), val_vectors, model, '-q');
    confusion{s} = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            confusion{s}(i,j) = nnz(answers(val_labels{s} == i) == j);
        end
    end
    cvacc(s) = sum(diag(confusion{s}))/sum(sum(confusion{s}));
end

fprintf('Average CV accuracy: %g\n', mean(cvacc));

%% plot results

% ridiculous reshaping to average the confusion matrices across CV splits
confall = cell2mat(cellfun(@(c) reshape(c, [69 1 69]), confusion, 'uniformoutput',false));
confavg = reshape(sum(confall, 2)/10, [69 69]);

% draw figure
% cut off "Gx" prefix for figure labels and remove tick marks
fig_confusion(confavg, cellfun(@(s) s(3:end), materials, 'uniformoutput',false), 5, 'Courier', 90, 0, false);
set(gca, 'TickLength', [0 0]);
print -dpdf cvconf.pdf;

