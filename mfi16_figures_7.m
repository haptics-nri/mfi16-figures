% part 7 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% setup for learning

% dataset parameters
materials = {'black', 'white', 'blue', 'brown', 'red'};
material_names = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
tools = {'stick'};
reps = {'1', '2', '3', '4', '5', '1', '2', '3', '4', '5'};
date = [repmat({'20160310'}, 1, 5) repmat({'20160219'}, 1, 5)];

%                 reps
% 2016-03-10 2016-02-19
video_offsets = [20  0 13 10  9 24 44 23 32 25
                 17 34 25 27 27 42 22 36 44 36 % materials
                 23 25 18 21 24 37 29 50 58 34
                 26 20 22 21 31 26 29 50 53 48
                 26 19 27 22 29 61 46 37 32 46];

% end-effector properties
mass = 0.1503; % kg
com = [-.0029; -.00301; .0348]; % m
