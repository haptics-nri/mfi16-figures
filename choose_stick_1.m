% part 1 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% weigh all sticks

datadir = '/Volumes/shared/Projects/Proton Pack/Data';
eps = {'20160620/weigh/1', '20160620/weigh/2', '20160617/weigh/1', '20160620/weigh/3'}; % TODO just scan the dir
sizes = [.25 .375 .5 .75];

masses = zeros(length(sizes),1);
coms = zeros(length(sizes),3);
bias = zeros(1,6);
for i=1:length(sizes)
    choose_stick_2
end
bias = bias / length(sizes);

endeff_idx = containers.Map;
endeff_idx('1/4') = 1;
endeff_idx('3/8') = 2;
endeff_idx('1/2') = 3;
endeff_idx('3/4') = 4;
material_idx = containers.Map;
material_idx('abs') = 1;
material_idx('plate') = 2;
material_idx('folder') = 3;
material_idx('mdf') = 4;
material_idx('canvas') = 5;
material_idx('weave') = 6;

save stickweights masses coms bias endeff_idx material_idx

