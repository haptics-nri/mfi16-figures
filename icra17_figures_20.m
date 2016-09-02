% part 20 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
%
% crossval

cv = cvpartition(cell2mat(train_features(:,1)), 'KFold', 3);
oc_confusion = cell(1, cv.NumTestSets);
mc_confusion = cell(1, cv.NumTestSets);
oc_answers = cell(1, cv.NumTestSets);
mc_answers = cell(1, cv.NumTestSets);
