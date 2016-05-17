% part 46 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% test on test set (first run GS/CV to completion)

[~, gsi] = max(gs_acc);
gs_nbins = nbins(gs_idx(gsi,1));
gs_binmode = binmode{gs_idx(gsi,2)};
gs_alpha = alpha(gs_idx(gsi,3));
gs_nu = nu(gs_idx(gsi,4));
gs_gamma = gamma(gs_idx(gsi,5));
gs_stmode = stmode(gs_idx(gsi,6));

train_vectors = [cell2mat(train_features(:,1)) ...
                          romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode)];
test_vectors  = [cell2mat(test_features(:,1)) ...
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
mc_train_args = [sprintf('-m 1000 -s 1 -t 2 -n %g -g %g', gs_nu, gs_gamma) common_args];
final_model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), mc_train_args);
mc_test_answers = svmpredict(zeros(size(test_vectors,1),1), test_vectors(:,2:end), final_model, test_args);
mc_test_confusion = zeros(length(materials));
for i=1:length(materials)
    mfi16_figures_47
end
fprintf('Test set accuracy: %g\n', sum(diag(mc_test_confusion))/sum(sum(mc_test_confusion)));

