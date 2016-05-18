% part 65 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
    fprintf('%s  &  ', material_names{i});
    others = setdiff(1:5, i);
    tp = mc_test_confusion(i,i);
    fp = sum(mc_test_confusion(others,i));
    tn = sum(sum(mc_test_confusion(others,others)));
    fn = sum(mc_test_confusion(i,others));
    acc = (tp+tn)/(tp+tn+fp+fn);
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    f1 = 2*prec*rec/(prec+rec);
    fprintf('%.3f  &  %.3f  &  %.3f  &  %.3f', acc, prec, rec, f1);
    fprintf(' \\\\ \n');
