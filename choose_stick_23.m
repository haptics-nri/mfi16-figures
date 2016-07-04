% part 23 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% cross-validation

% hyperparameters
nbins = 1:2:20; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.5; % 25
nu = .05:0.05:0.7; % .6
gamma = 1:2:20; % 200
stmode = true; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for j=2:size(gs_idx,1)
    choose_stick_24
end

for i=1:4
    choose_stick_28
end

save endeffdata endeffs

