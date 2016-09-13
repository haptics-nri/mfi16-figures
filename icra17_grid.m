function gs_acc = icra17_grid(materials, features, split_idx, nbins, binmode, alpha, nu, gamma)

    train_features = features(split_idx==1, :);
    
    cv = cvpartition(cell2mat(train_features(:,1)), 'KFold', 3);
    mc_confusion = cell(1, cv.NumTestSets);
    mc_answers = cell(1, cv.NumTestSets);
    
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
    
    gs_acc = zeros(size(gs_idx,1),1);
    clear romano_features; % clear persistent vars
    elapsed = tic;
    for gsi=1:size(gs_idx,1)
        gs_nbins = nbins(gs_idx(gsi,1));
        gs_binmode = binmode{gs_idx(gsi,2)};
        gs_alpha = alpha(gs_idx(gsi,3));
        gs_nu = nu(gs_idx(gsi,4));
        gs_gamma = gamma(gs_idx(gsi,5));
        
        iter = tic;
        fprintf('Grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g\n', gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma);
        cv_acc = zeros(cv.NumTestSets,1);
        for cvi = 1:cv.NumTestSets
            train_vectors = [cell2mat(train_features(cv.training(cvi),1)) ...
                             romano_features('post', train_features(cv.training(cvi),2:end), gs_nbins, gs_binmode, gs_alpha, 1)];
            val_vectors   = [cell2mat(train_features(cv.test(cvi),    1)) ...
                             romano_features('post', train_features(cv.test(cvi)    ,2:end), gs_nbins, gs_binmode, gs_alpha, 1)];

            trainmean = mean(train_vectors(:,2:end));
            train_vectors(:,2:end) = bsxfun(@minus, ...
                                            train_vectors(:,2:end), ...
                                            trainmean);
            val_vectors  (:,2:end) = bsxfun(@minus, ...
                                            val_vectors  (:,2:end), ...
                                            trainmean);
            trainrange = range(train_vectors(:,2:end));
            train_vectors(:,2:end) = bsxfun(@rdivide, ...
                                            train_vectors(:,2:end), ...
                                            trainrange);
            val_vectors  (:,2:end) = bsxfun(@rdivide, ...
                                            val_vectors  (:,2:end), ...
                                            trainrange);

            mc_train_args = sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', gs_nu, gs_gamma);

            model = svmtrain(train_vectors(:,1), train_vectors(:,2:end), mc_train_args);

            if ~isempty(model)
                mc_confusion{cvi} = zeros(length(materials));
                mc_answers{cvi} = svmpredict(zeros(size(val_vectors,1),1), val_vectors(:,2:end), model, '-q');

                for i=1:length(materials)
                    for j=1:length(materials)
                        mc_confusion{cvi}(i,j) = nnz(mc_answers{cvi}(val_vectors(:,1)==i) == j);
                    end
                end

                cv_acc(cvi) = sum(diag(mc_confusion{cvi}))/sum(sum(mc_confusion{cvi}));
                fprintf('\tFold %d: MC %g%%\n', cvi, 100*cv_acc(cvi));
            else
                fprintf('\tFold %d: failed to train\n', cvi);
                cv_acc(cvi) = 0;
            end
        end
        gs_acc(gsi) = mean(cv_acc);
        fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi), toc(iter), toc(elapsed));
    end

end
