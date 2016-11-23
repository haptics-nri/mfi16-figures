% part 8 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% setup for learning

% dataset parameters
date = '20160906';
flowtype = 'stickcam';

% end-effector mass comes from calibration above

[data14, materials14] = icra17_load(DATADIR, date, flowtype, @(x) x < 7);
[data38, materials38] = icra17_load(DATADIR, date, flowtype, @(x) x > 5);
assert(all(cellfun(@strcmp, materials14, materials38)));
materials = materials14;

