% part 17 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    d{i}.pre = romano_features('pre', d{i}.iws, d{i}.vei, d{i}.ai, mass, 150, 0, [d{i}.a-1515 d{i}.b+1515]);
    d{i}.pre = d{i}.pre(11:20,:);
    d{i}.post = romano_features('post', d{i}.pre, 5, 'naive', 0, 1);
    d{i}.post(:, (end-3):end) = d{i}.post(:,[end-1 end end-3 end-2]); % swap V and Ft
    
    d{i}.m = mean(d{i}.post);
    d{i}.r = range(bsxfun(@minus, d{i}.post, d{i}.m));
    d{i}.m(1:5) = mean(d{i}.m(1:5));
    d{i}.r(1:5) = mean(d{i}.r(1:5))*3;
    d{i}.m([6 8]) = mean(d{i}.m([6 8]));
    d{i}.r([6 8]) = mean(d{i}.r([6 8]));
