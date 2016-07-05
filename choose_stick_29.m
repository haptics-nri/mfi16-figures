% part 29 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        gs_nbins = nbins(gs_idx(gsi,1));
        gs_binmode = binmode{gs_idx(gsi,2)};
        gs_alpha = alpha(gs_idx(gsi,3));
        gs_nu = nu(gs_idx(gsi,4));
        gs_gamma = gamma(gs_idx(gsi,5));
        gs_stmode = stmode(gs_idx(gsi,6));
        
        iter = tic;
        fprintf('endeff %d grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g, stmode=%d\n', i, gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma, gs_stmode);
        cv_acc = zeros(cv.NumTestSets,1);
        for cvi = 1:cv.NumTestSets
            choose_stick_30
        end
        gs_acc(gsi) = mean(cv_acc);
        endeffs(i).gs.cv_acc{gsi} = cv_acc;
        fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi), toc(iter), toc(elapsed));
