% part 31 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
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
                choose_stick_32
            else
                choose_stick_35
            end
