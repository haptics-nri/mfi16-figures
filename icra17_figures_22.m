% part 22 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
%
% hyperparameters
nbins = 10:10:60; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.5; % 25
nu = .05:0.05:0.3; % .6
gamma = 10:20:100; % 200
stmode = true; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for i=2:size(gs_idx,1)
    icra17_figures_23
end
