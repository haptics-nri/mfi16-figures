% part 24 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    train_features = endeffs(i).features(endeffs(i).split==1, :);
    test_features  = endeffs(i).features(endeffs(i).split==2, :);

    train_vectors = [cell2mat(train_features(:,1)) ...
                     romano_features('post', train_features(:,2:end), 3, 'perceptual', 0.2, true)];
    test_vectors  = [cell2mat(test_features (:,1)) ...
                     romano_features('post', test_features (:,2:end), 3, 'perceptual', 0.2, true)];

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

    common_args = ' -q ';
    train_args = [sprintf('-m 1000 -s 1 -t 2 -n %g -g %g', 0.1, 7) common_args];
    test_args = common_args;

    model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), train_args);

    if ~isempty(model)
        choose_stick_25
    else
        choose_stick_28
    end
