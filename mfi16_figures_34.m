% part 34 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
%
% hyperparameters
nbins = 5:5:40; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.05:0.05:0.5; % 25
nu = .2:0.05:0.6; % .6
gamma = [1 5 10 20]; % 200
stmode = [false true]; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for i=2:size(gs_idx,1)
    mfi16_figures_35
end
