% part 3 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    cc = xfconv(x(i,2:7)) \ [c 1]';
    d(i,:) = cc(1:3);
