% part 44 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    [~,gsi] = max(endeffs(i).gs_acc);
    gs_nbins = nbins(gs_idx(gsi,1));
    gs_binmode = binmode{gs_idx(gsi,2)};
    gs_alpha = alpha(gs_idx(gsi,3));
    gs_nu = nu(gs_idx(gsi,4));
    gs_gamma = gamma(gs_idx(gsi,5));
    gs_stmode = stmode(gs_idx(gsi,6));
    
    train_features = endeffs(i).features(endeffs(i).split==1, :);
    test_features = endeffs(i).features(endeffs(i).split==2, :);
            
    train_vectors = [cell2mat(train_features(:,1)) ...
                     romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];
    test_vectors  = [cell2mat(test_features (:,1)) ...
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
    
    endeffs(i).model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), sprintf('-q -m 1000 -s 1 -t 2 -n %g -g %g', gs_nu, gs_gamma));
    endeffs(i).predictions = svmpredict(zeros(size(test_vectors,1),1), test_vectors(:,2:end), endeffs(i).model, '-q');
    endeffs(i).confusion = zeros(5,5);
    for j=1:5
        choose_stick_45
    end
    endeffs(i).accuracy = sum(diag(endeffs(i).confusion))/sum(sum(endeffs(i).confusion));
