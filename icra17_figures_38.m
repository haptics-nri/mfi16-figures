% part 38 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% test on test set (first run GS/CV to completion)

[~, gsi] = max(gs_acc);
gs_nbins = nbins(gs_idx(gsi,1));
gs_binmode = binmode{gs_idx(gsi,2)};
gs_alpha = alpha(gs_idx(gsi,3));
gs_nu = nu(gs_idx(gsi,4));
gs_gamma = gamma(gs_idx(gsi,5));
gs_stmode = stmode(gs_idx(gsi,6));
