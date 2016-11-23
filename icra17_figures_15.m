% part 15 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% feature vector figure

dosub = false;
if dosub
    icra17_figures_16
end

% sample data
d1 = data38('abs');
d1.vel = pose_to_vel(d1.vei, d1.iws);
d1.a = 47200;
d1.b = d1.a+1500;
d{1} = d1;
d2 = data38('vinyl');
d2.vel = pose_to_vel(d2.vei, d2.iws);
d2.a = 55200;
d2.b = d2.a+1500;
d{2} = d2;

for i=1:2
    icra17_figures_17
end

meanm = mean([d{1}.m; d{2}.m]);
meanr = mean([d{1}.r; d{2}.r]);
for i=1:2
    icra17_figures_18
end
minmin = min(minmin);
maxmax = max(maxmax);

f = figure;
imagesc([d{1}.posted d{2}.posted]);
ax = gca;
clim = ax.CLim;
close(f);

for i=1:2
    icra17_figures_19
end

if ~dosub
    icra17_figures_31
end
