% part 60 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% GS sensitivity analysis

gs_plots(vgs_acc, gs_idx(:,[1 3 4 5]), {'nbins', nbins; 'alpha', alpha; 'nu', nu; 'gamma', gamma});
suplabel('Vicon', 't');

gs_plots(bgs_acc, gs_idx(:,[1 3 4 5]), {'nbins', nbins; 'alpha', alpha; 'nu', nu; 'gamma', gamma});
suplabel('Bluefox', 't');
