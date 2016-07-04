%% weigh all sticks

datadir = '/Volumes/shared/Projects/Proton Pack/Data';
eps = {'20160620/weigh/1', '20160620/weigh/2', '20160617/weigh/1', '20160620/weigh/3'}; % TODO just scan the dir
sizes = [.25 .375 .5 .75];

masses = zeros(length(sizes),1);
coms = zeros(length(sizes),3);
bias = zeros(1,6);
for i=1:length(sizes)
    [masses(i), fbias, ~, coms(i,:), tbias, ~] = weigh({[datadir filesep eps{i}]});
    bias(1:3) = bias(1:3) + fbias;
    bias(4:6) = bias(4:6) + tbias';
end
bias = bias / length(sizes);

endeff_idx = containers.Map;
endeff_idx('1/4') = 1;
endeff_idx('3/8') = 2;
endeff_idx('1/2') = 3;
endeff_idx('3/4') = 4;
material_idx = containers.Map;
material_idx('abs') = 1;
material_idx('plate') = 2;
material_idx('folder') = 3;
material_idx('mdf') = 4;
material_idx('canvas') = 5;
material_idx('weave') = 6;

save stickweights masses coms bias endeff_idx material_idx

%% load texture data

date = '20160620';
epdirs = dir([datadir filesep date filesep 'stick']);
epdirs(arrayfun(@(e) e.name(1) == '.', epdirs)) = [];

eps = struct('endeff', cell(length(epdirs),1), ...
             'material', cell(length(epdirs),1), ...
             'flow', cell(length(epdirs),1), ...
             'data', cell(length(epdirs),1), ...
             'features', cell(length(epdirs),1));
for i=1:length(epdirs)
    fprintf('[%d/%d] %s\n', i, length(epdirs), epdirs(i).name);
    j = str2double(epdirs(i).name);
    prefix = [datadir filesep date filesep 'stick' filesep epdirs(i).name];
    eps(j).flow = parse_flow([prefix filesep 'stick.flow']);
    eps(j).endeff = eps(j).flow.answers('tooling ball diameter').text;
    eps(j).material = eps(j).flow.answers('surface name').text;
    [v, f, ~,~,~, a] = load_stick([prefix filesep]);
    eps(j).data = struct('vicon', v, 'force', f, 'acc', a);
end

save stickdata eps

%% sync vicon/force data

for i=1:length(eps)
    fprintf('%d\n', i);
    
    v = eps(i).data.vicon;
    f = eps(i).data.force;
    
    % find the 4 taps in force (hand-tuned parameters seem to work)
    [fpks, flocs] = findpeaks(f(:,3), 'SortStr','descend', 'NPeaks',4, 'MinPeakProminence',2, 'MinPeakDistance',500);
    [flocs, idx] = sort(flocs); % resort in order
    fpks = fpks(idx);
    ftimes = f(flocs([1 3]),1); % use the first and second-to-last
    
    % find the 4 rises-before-the-taps in vicon
    [vpks, vlocs] = findpeaks(v(:,4), 'SortStr','descend', 'NPeaks',4, 'MinPeakDistance',10);
    [vlocs, idx] = sort(vlocs); % resort in order
    vpks = vpks(idx);
    [~,i1] = min(v(vlocs(1):vlocs(2),4));
    [~,i2] = min(v(vlocs(3):vlocs(4),4));
    vtimes = [v(i1 + vlocs(1) - 1, 1); v(i2 + vlocs(3) - 1, 1)]; % find min between each pair to get the first and second-to-last
    
    % the offset is the average time difference between the peaks,
    % but we might have some misidentified peaks.
    eps(i).offset = mean(vtimes - ftimes);
    eps(i).peaks = flocs;
end
% loop through and fixup any that are too far from the mean
for idx = find(abs([eps.offset] - mean([eps.offset])) > 0.25)
    eps(idx).offset = (eps(idx-1).offset + eps(idx+1).offset)/2;
end

save stickdata eps

%% process data

H_vic2bod = ...
   [0.9912   -0.0236   -0.1303         0
    0.0162    0.9982   -0.0571   36.0685
    0.1314    0.0545    0.9898 -511.6330
         0         0         0    1.0000];
    
 H_bal2imu = ...
   [1.0000         0         0  254.3402
         0    1.0000         0         0
         0         0    1.0000         0
         0         0         0    1.0000];
     
 H_m402bod = ...
   [     0         0    1.0000  108.9900
    1.0000         0         0    0.5300
         0    1.0000         0   -2.9800
         0         0         0    1.0000];

for i=1:length(eps)
    fprintf('%d\n', i);
    if std(diff(eps(i).data.force(100:end,1))) < 2e-4
        if ~isfield(eps(i).data, 'f_comp')
            j = endeff_idx(eps(i).endeff);
            [~,~,~,~,~,~, eps(i).data.v_end, eps(i).data.a_end, ~,~,~, eps(i).data.f_comp] = process_stick(eps(i).data.vicon, eps(i).data.force, eps(i).data.acc, masses(j), coms(j,:)', H_vic2bod, H_m402bod, H_bal2imu, -eps(i).offset);
        end
    end
end

save stickdata eps

%% some plots

figs = zeros(1,4);
order = 1:6;
for e=1:4
    figs(e) = figure;
    for s=1:length(order)
        subplot(2,3,s);
        k = (e-1)*30 + (order(s)-1)*5 + 1;
        if ~isfield(eps(k).data, 'f_comp')
            k = k + 1;
        end
        a = eps(k).peaks(2) + 1500;
        b = eps(k).peaks(3) - 1500;
        plot(eps(k).data.f_comp(a:b,1)-eps(k).data.f_comp(1,1), eps(k).data.f_comp(a:b,2:4), ...
             eps(k).data.v_end(a:b,1)-eps(k).data.f_comp(1,1), bsxfun(@minus, eps(k).data.v_end(a:b,2:3), eps(k).data.v_end(a,2:3))/10 + 20);
        title(eps(k).material);
        axis tight;
    end
    suplabel(eps(k).endeff, 't');
    print('-depsc', sprintf('choose_stick__%d_eighths.eps', eval(eps(k).endeff)*8));
    order = fliplr(order);
end

%% run classification with same parameters as in the paper

for i=1:length(eps)
    fprintf('%d calculating pre-features... ', i);
    j = endeff_idx(eps(i).endeff);
    if isfield(eps(i).data, 'f_comp')
        eps(i).features.pre = romano_features('pre', eps(i).data.f_comp, eps(i).data.v_end, eps(i).data.a_end, masses(j), 0.05, [20 3], [eps(i).peaks(2)+1000 eps(i).peaks(3)-1000]);
        eps(i).features.post = romano_features('post', eps(i).features.pre, 3, 'perceptual', 0.2, true);
    end
    fprintf('done\n');
end

%save stickdata eps

%% test/train split

endeffs = struct('features', cell(4,1), ...
                 'split', cell(4,1), ...
                 'gs', cell(4,1), ...
                 'model', cell(4,1), ...
                 'predictions', cell(4,1), ...
                 'confusion', cell(4,1), ...
                 'accuracy', cell(4,1));

for i=1:length(eps)
    j = endeff_idx(eps(i).endeff);
    if isfield(eps(i).data, 'f_comp') && ~strcmp(eps(i).material, 'weave')
        endeffs(j).features = [endeffs(j).features
                               num2cell(repmat(material_idx(eps(i).material), [size(eps(i).features.pre,1) 1])) eps(i).features.pre];
    end
end
   
for i=1:4
    % 4/5 train, 1/5 test
    endeffs(i).split = randsample(1:2, size(endeffs(i).features,1), true, [4/5 1/5]);
end

save endeffdata endeffs

%% cross-validation

% hyperparameters
nbins = 1:2:20; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.5; % 25
nu = .05:0.05:0.7; % .6
gamma = 1:2:20; % 200
stmode = true; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for j=2:size(gs_idx,1)
    gs_idx(j,:) = gs_idx(j-1,:);
    for k=size(gs_idx,2):-1:1
        if gs_idx(j,k) == gs_limits(k)
            gs_idx(j,k) = 1;
        else
            gs_idx(j,k) = gs_idx(j,k) + 1;
            break;
        end
    end
end

for i=1:2
    % CV partition
    cv = cvpartition(cell2mat(endeffs(i).features(endeffs(i).split==1, 1)), 'KFold', 3);
    confusion = cell(1, cv.NumTestSets);
    predictions = cell(1, cv.NumTestSets);
    endeffs(i).gs.partition = cv;
    endeffs(i).gs.cv_acc = cell(size(gs_idx,1),1);
    endeffs(i).gs.gs_acc = zeros(size(gs_idx,1),1);
    
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
        
        iter = tic;
        fprintf('endeff %d grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g, stmode=%d\n', i, gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma, gs_stmode);
        cv_acc = zeros(cv.NumTestSets,1);
        for cvi = 1:cv.NumTestSets
            train_features = endeffs(i).features(endeffs(i).split==1, :);
            
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

            common_args = ' -q ';
            train_args = [sprintf('-m 1000 -s 1 -t 2 -n %g -g %g', gs_nu, gs_gamma) common_args];
            test_args = common_args;

            model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), train_args);

            if ~isempty(model)
                confusion{cvi} = zeros(5,5);
                predictions{cvi} = svmpredict(zeros(size(val_vectors,1),1), val_vectors(:,2:end), model, '-q');

                for j=1:5
                    for k=1:5
                        confusion{cvi}(j,k) = nnz(predictions{cvi}(val_vectors(:,1)==j) == k);
                    end
                end

                cv_acc(cvi) = sum(diag(confusion{cvi}))/sum(sum(confusion{cvi}));
                fprintf('\tFold %d: MC %g%%\n', cvi, 100*cv_acc(cvi));
            else
                fprintf('\tFold %d: failed to train\n', cvi);
                cv_acc(cvi) = 0;
            end
        end
        gs_acc(gsi) = mean(cv_acc);
        endeffs(i).gs.cv_acc{gsi} = cv_acc;
        fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi), toc(iter), toc(elapsed));
    end
    endeffs(i).gs_acc = gs_acc;
end

save endeffdata endeffs

%% run classification

for i=1:4
    endeffs(i).model = svmtrain(endeffs(i).vectors.train(:,1), endeffs(i).vectors.train(:,2:end), '-q -m 1000 -s 1 -t 2 -n 0.1 -g 7.0');
    endeffs(i).predictions = svmpredict(zeros(size(endeffs(i).vectors.test,1),1), endeffs(i).vectors.test(:,2:end), endeffs(i).model, '-q');
    endeffs(i).confusion = zeros(5,5);
    for j=1:5
        for k=1:5
            endeffs(i).confusion(j,k) = nnz(endeffs(i).predictions(endeffs(i).vectors.test(:,1)==j) == k);
        end
    end
    endeffs(i).accuracy = sum(diag(endeffs(i).confusion))/sum(sum(endeffs(i).confusion));
end

%save endeffdata endeffs

