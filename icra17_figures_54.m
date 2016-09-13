% part 54 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    fprintf('%s  &  ', material_names{i});
    others = setdiff(1:5, i);
    tp = mc_test_confusion(i,i);
    fp = sum(mc_test_confusion(others,i));
    tn = sum(sum(mc_test_confusion(others,others)));
    fn = sum(mc_test_confusion(i,others));
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    f1 = 2*prec*rec/(prec+rec);
    fprintf('%.3f  &  %.3f  &  %.3f', prec, rec, f1);
    fprintf(' \\\\ \n');
