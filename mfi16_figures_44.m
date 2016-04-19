% part 44 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
%
clf;
subplot(1,2,1);
bar3(mc_confusion);
title(sprintf('MC accuracy = %g%%', 100*sum(diag(mc_confusion))/sum(sum(mc_confusion))));
subplot(1,2,2);
bar3(oc_confusion);
title(sprintf('OC accuracy = %g%%', 100*sum(diag(oc_confusion))/sum(sum(oc_confusion))));

