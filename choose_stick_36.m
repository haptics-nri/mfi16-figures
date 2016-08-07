% part 36 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
       %%
        iter = tic;
        fprintf('endeff %d grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g, stmode=%d\n', i, gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma, gs_stmode);
        cv_acc = zeros(cv.NumTestSets,1);
        for cvi = 1:cv.NumTestSets
            choose_stick_37
        end
        fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*mean(cv_acc), toc(iter), toc(elapsed));
