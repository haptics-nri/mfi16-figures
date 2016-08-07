% part 29 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% cross-validation

% hyperparameters
nbins = 7:2:20; % 20
binmode = {'perceptual'}; % perceptual
alpha = 0.1:0.05:0.4; % 25
nu = .1:0.05:0.3; % .6
gamma = 1:.5:5; % 200
stmode = true; % false

gs_limits = [length(nbins) length(binmode) length(alpha) length(nu) length(gamma) length(stmode)];
gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
for j=2:size(gs_idx,1)
    choose_stick_30
end

for i=1:2
    choose_stick_34
end

%save endeffdata endeffs

