% Grid search for WHC 17
function gs_acc = whc17_grid(materials, labels, features, split_idx, nbins, binmode, alpha, nu, gamma)

    % pull out training set
    labels = labels(split_idx==1);
    features{1} = features{1}(split_idx==1, :);
    features{2} = features{2}(split_idx==1, :);
    
    % prepare cross validation
    cv = cvpartition(labels, 'KFold', 3);
    
    % prepare grid search
    gs_idx = prepare_grid(nbins, binmode, alpha, nu, gamma);
    
    % execute grid search!
    gs_acc = zeros(size(gs_idx,1),3);
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
        cv_acc = zeros(cv.NumTestSets,3);
        for cvi = 1:cv.NumTestSets
            train_labels = labels(cv.training(cvi));
            val_labels   = labels(cv.test(cvi)    );
            train_vectors = [romano_features('post', features{1}(cv.training(cvi),:), gs_nbins, gs_binmode, gs_alpha, 1) ...
                             features{2}(cv.training(cvi),:)];
            val_vectors   = [romano_features('post', features{1}(cv.test(cvi)    ,:), gs_nbins, gs_binmode, gs_alpha, 1) ...
                             features{2}(cv.test(cvi)    ,:)];

            trainmean = mean(train_vectors);
            train_vectors = bsxfun(@minus, ...
                                   train_vectors, ...
                                   trainmean);
            val_vectors   = bsxfun(@minus, ...
                                   val_vectors  , ...
                                   trainmean);
            trainrange = range(train_vectors);
            train_vectors = bsxfun(@rdivide, ...
                                   train_vectors, ...
                                   trainrange);
            val_vectors   = bsxfun(@rdivide, ...
                                   val_vectors  , ...
                                   trainrange);

            mc_train_args = sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', gs_nu, gs_gamma);

            % Romano-only, Steinbach-only, and both
            slen = size(features{2}(cv.training(cvi),:), 2);
            cv_acc(cvi,1) = do_cv(train_labels, train_vectors(:,1:end-slen), val_labels, val_vectors(:,1:end-slen), materials, mc_train_args);
            cv_acc(cvi,2) = do_cv(train_labels, train_vectors(:,end-slen+1:end), val_labels, val_vectors(:,end-slen+1:end), materials, mc_train_args);
            cv_acc(cvi,3) = do_cv(train_labels, train_vectors, val_labels, val_vectors, materials, mc_train_args);
            fprintf('\tFold %d: r=%g, s=%g, rs=%g\n', cvi, cv_acc(cvi,1), cv_acc(cvi,2), cv_acc(cvi,3));
        end
        gs_acc(gsi,:) = mean(cv_acc);
        fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi,3), toc(iter), toc(elapsed));
    end

end

function cv_acc = do_cv(train_labels, train_vectors, val_labels, val_vectors, materials, mc_train_args)
    model = svmtrain(train_labels, train_vectors, mc_train_args);

    if ~isempty(model)
        mc_confusion = zeros(length(materials));
        mc_answers = svmpredict(zeros(size(val_vectors,1),1), val_vectors, model, '-q');

        for i=1:length(materials)
            for j=1:length(materials)
                mc_confusion(i,j) = nnz(mc_answers(val_labels==i) == j);
            end
        end

        cv_acc = sum(diag(mc_confusion))/sum(sum(mc_confusion));
    else
        cv_acc = 0;
    end
end
