function icra17_test(materials, features, split_idx, gs_acc, nbins, binmode, alpha, nu, gamma)

    train_features = features(split_idx==1, :);
    test_features = features(split_idx==2, :);
    
    gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma)];
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
    [~, gsi] = max(gs_acc);
    nbins = nbins(gs_idx(gsi,1));
    binmode = binmode{gs_idx(gsi,2)};
    alpha = alpha(gs_idx(gsi,3));
    nu = nu(gs_idx(gsi,4));
    gamma = gamma(gs_idx(gsi,5));
    
    fprintf('Test set with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g (params from max GS=%g)\n', nbins, binmode, alpha, nu, gamma, max(gs_acc));

    train_vectors = [cell2mat(train_features(:,1)) ...
                              romano_features('post', train_features(:,2:end), nbins, binmode, alpha, 1)];
    test_vectors  = [cell2mat(test_features(:,1)) ...
                              romano_features('post', test_features (:,2:end), nbins, binmode, alpha, 1)];
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
    mc_train_args = sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', nu, gamma);
    final_model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), mc_train_args);
    mc_test_answers = svmpredict(zeros(size(test_vectors,1),1), test_vectors(:,2:end), final_model, '-q');
    mc_test_confusion = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            mc_test_confusion(i,j) = nnz(mc_test_answers(test_vectors(:,1)==i) == j);
        end
    end
    fprintf('Test set accuracy: %g\n', sum(diag(mc_test_confusion))/sum(sum(mc_test_confusion)));

end