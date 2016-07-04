% part 2 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    [masses(i), fbias, ~, coms(i,:), tbias, ~] = weigh({[datadir filesep eps{i}]});
    bias(1:3) = bias(1:3) + fbias;
    bias(4:6) = bias(4:6) + tbias';
