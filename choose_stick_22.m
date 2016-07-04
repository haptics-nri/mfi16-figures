% part 22 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    % 4/5 train, 1/5 test
    endeffs(i).split = randsample(1:2, size(endeffs(i).features,1), true, [4/5 1/5]);
