% part 11 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% SVM stuff (following Romano & KJK 2014 + Strese & Schuwerk & Steinbach 2015)

[features14, split_idx14, bfeatures14, bsplit_idx14] = icra17_svm(data14, mass);
[features38, split_idx38, bfeatures38, bsplit_idx38] = icra17_svm(data38, mass);

