% part 42 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    idx = cell2mat(train_features(:,1))==i;
    g = f(idx,:);
    g = bsxfun(@minus, g, mean(g));
    g = bsxfun(@rdivide, g, range(g));
    g = [g mean(g(:,[end-5 end-3 end-1]), 2)];
    g = sortrows(g, size(g,2));
    g = g(:,1:end-1);
    f(idx,:) = g;
