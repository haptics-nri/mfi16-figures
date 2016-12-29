%% preprocess

for m=1:length(materials)
    fprintf('%s\n', materials{m});
    ep = data14(materials{m});
    n = floor(diff(ep.ss)/3000);
    ep.sf = zeros(n,16);
    for i=1:n
        ep.sf(i,:) = steinbach_features({'MF', 'SR', 'F', 'RG'}, [0 0], [ep.ss(1)+(i-1)*3000 ep.ss(1)+i*3000], ...
                                        struct('data', dft321(ep.acc(:,[2 4 5])), 'Fs', 3000), [], []);
    end
    data14(materials{m}) = ep;
    
    ep = data14new(materials{m});
    n = floor(size(ep.acc,1)/3000);
    ep.sf = zeros(n,16);
    for i=1:n
        ep.sf(i,:) = steinbach_features({'MF', 'SR', 'F', 'RG'}, [0 0], [1 size(ep.acc,1)], ...
                                        struct('data', dft321(ep.acc(:,[2 4 5])), 'Fs', 3000), [], []);
    end
    data14new(materials{m}) = ep;
end

sf = cell2mat(cellfun(@(s, i) [repmat(i,size(s.sf,1),1) s.sf], data14.values, num2cell(1:length(materials)), 'uniformoutput',false)');
sfnew = cell2mat(cellfun(@(s, i) [repmat(i,size(s.sf,1),1) s.sf], data14new.values, num2cell(1:length(materials)), 'uniformoutput',false)');
sf = [sf; sfnew];

%% divide

train_features = [];
val_features = [];
train_labels = [];
val_labels = [];

for i=1:length(materials)
    idx = find(sf(:,1)==i);
    d = randsample(idx, round(.5*length(idx)));
    train_labels = [train_labels; repmat(i, length(d),1)];
    train_features = [train_features; sf(d,2:end)];
    d = setdiff(idx, d);
    val_labels = [val_labels; repmat(i, length(d),1)];
    val_features = [val_features; sf(d,2:end)];
end

%% data from separate days

train_labels = sf(:,1);
train_features = sf(:,2:end);
val_labels = sfnew(:,1);
val_features = sfnew(:,2:end);

%% test

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

model = svmtrain(train_labels, train_vectors, '-m 1000 -s 1 -t 2 -n 0.1 -g 0.1 -q');
answers = svmpredict(zeros(size(val_labels)), val_vectors, model, '-q');
confusion = zeros(length(materials));
for i=1:length(materials)
    for j=1:length(materials)
        confusion(i,j) = nnz(answers(val_labels == i) == j);
    end
end
fprintf('Test set accuracy: %g\n', sum(diag(confusion))/sum(sum(confusion)));
