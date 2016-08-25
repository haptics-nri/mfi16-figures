% part 8 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% setup for learning

% dataset parameters
date = '20160815';
flowtype = 'stickcam';

% end-effector mass comes from calibration above

data = containers.Map;
episodes = dir([DATADIR filesep date filesep flowtype]);
for ep = 1:length(episodes)
    icra17_figures_9
end

materials = data.keys;

for m = 1:length(materials)
    icra17_figures_11
end

