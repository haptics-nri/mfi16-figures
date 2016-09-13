% part 8 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% setup for learning

% dataset parameters
date = '20160906';
flowtype = 'stickcam';

% end-effector mass comes from calibration above

data = containers.Map;
episodes = dir([DATADIR filesep date filesep flowtype]);
for ep = 1:length(episodes)
    icra17_figures_9
end

materials = data.keys;

if false
    icra17_figures_12
end

