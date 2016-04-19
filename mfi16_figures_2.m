% part 2 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% sphere calibration (see go_sphere_again.m)

% vicon data
v2 = csvload([DATADIR '/20160223/socket2stick/vicon.tsv'], ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});
v3 = csvload([DATADIR '/20160223/socket3stick/vicon.tsv'], ...
             {'Timestamp', ...
              'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
              'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
             {'Delimiter', '\t'});

% use the second and third (ball popped up in first)
x = [v2; v3];

[c, r, ~, cs] = sphereFit_ransac(x(:,2:4), 50); % FIXME allowing a huge amount of noise to get a reasonable number of inliers

d = nan([size(x,1) 3]);
for i = find(cs)
    mfi16_figures_3
end
d = nanmean(d);

% the product of sphere calibration is d
spherecalib.x = x;
spherecalib.c = c;
spherecalib.r = r;
spherecalib.cs = cs;
spherecalib.d = d;

