function [conf, acc, model] = go(materials, train_labels, train_features, test_labels, test_features, nu, gamma)
    trainmean = mean(train_features);
    trainrange = range(train_features);
    train_vectors = bsxfun(@rdivide, bsxfun(@minus, train_features, trainmean), trainrange);
    test_vectors = bsxfun(@rdivide, bsxfun(@minus, test_features, trainmean), trainrange);
    
    model = svmtrain(train_labels, train_vectors, sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', nu, gamma));
    answers = svmpredict(zeros(size(test_labels)), test_vectors, model, '-q');
    conf = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            conf(i,j) = nnz(answers(test_labels == i) == j);
        end
    end
    acc = sum(diag(conf))/sum(sum(conf));
end

