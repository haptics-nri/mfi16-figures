% part 29 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
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
            icra17_figures_30
        else
            icra17_figures_33
        end
