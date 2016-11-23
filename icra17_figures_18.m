% part 18 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    d{i}.posted = bsxfun(@rdivide, bsxfun(@minus, d{i}.post, meanm), meanr);
    minmin(i) = min(min(d{i}.posted));
    maxmax(i) = max(max(d{i}.posted));
