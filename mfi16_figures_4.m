% part 4 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% free-space calibration (see go_bias.m)

date = '20160226'; %'20151218';
offset = 22.848; % manual (from taps)
material = 'free';
tool = 'stick';
rep = '1';

set = [material rep tool];
prefix = [DATADIR filesep date filesep set filesep];

v = csvload([prefix 'vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

int = csvload([prefix 'teensy.ft.csv'], ...
              [{'Timestamp'}, ...
               arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
           
int = process_mini40(int, zeros(1,6));

a = round(size(int,1)*1/5);
b = round(size(int,1)*4/5);
int_s = int(a:b,:);

[~,aa] = min(abs(v(:,1)-int(a,1)));
[~,bb] = min(abs(v(:,1)-int(b,1)));
v_s = v(aa:bb,:);

[c,r,e,inliers] = sphereFit_ransac(int_s(:,2:4)); % 94% inliers

[~, start] = min(abs(int_s(:,1) - v(1,1)));

R_m402bod = [ 0       0       1
              1       0       0
              0       1       0 ];

grav = zeros(size(int_s,1), 4);
grav(:,1) = int_s(:,1) - offset;
ideal_grav = grav;
for tm=1:size(grav,1)
    mfi16_figures_5
end

[R, t, err] = rigid_ransac(ideal_grav(inliers,2:4), grav(inliers,2:4), 10); % 100% inliers??

% the product of free-space calibration is R
freecalib.int_s = int_s;
freecalib.grav = grav;
freecalib.ideal_grav = ideal_grav;
freecalib.c = c;
freecalib.r = r;
freecalib.cs = inliers;
freecalib.R = R;
freecalib.t = t;

