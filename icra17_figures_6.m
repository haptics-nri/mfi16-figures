% part 6 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% end-effector weighing

date = '20160811';
[~,f] = load_stick([DATADIR filesep date filesep 'weigh/1/']);
[~,r] = sphereFit_ransac(f(:,2:4)); % 98.9% inliers

mass = r/9.81;

