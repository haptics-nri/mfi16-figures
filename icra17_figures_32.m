% part 32 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
   %%
    iter = tic;
    fprintf('Grid search with nbins=%d, binmode=%s, alpha=%g, nu=%g, gamma=%g, stmode=%d\n', gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma, gs_stmode);
    cv_acc = zeros(cv.NumTestSets,1);
    for cvi = 1:cv.NumTestSets
        icra17_figures_33
    end
    gs_acc(gsi) = mean(cv_acc);
    fprintf('\tGS #%d/%d: mean acc %g%%, iter %g s / elapsed %g s\n', gsi, size(gs_idx,1), 100*gs_acc(gsi), toc(iter), toc(elapsed));
