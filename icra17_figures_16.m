% part 16 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% SVM stuff (following Romano & KJK 2014 + Strese & Schuwerk & Steinbach 2015)

% DO NOT RUN THIS AGAIN -- TRAINING ON TEST SET WILL RESULT

% extract features
features = cell(0, 5); % cols: label, vibration, speed, normal, tangential
for m = 1:length(materials)
    icra17_figures_17
    %%
    icra17_figures_18
    %%
    icra17_figures_19
end

% test/train split

% 4/5 train, 1/5 test
split_idx = randsample(1:2, size(features,1), true, [4/5 1/5]);

train_features = features(split_idx==1, :);
test_features  = features(split_idx==2, :);

