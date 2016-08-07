% part 3 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% load texture data

clear dir;
datadir = '/Volumes/shared/Projects/Proton Pack/Data';
date = '20160620';
epdirs = dir([datadir filesep date filesep 'stick']);
[epdirs.dir] = deal(date);
date = '20160713'; % we re-did the 3/4 trials
new_epdirs = dir([datadir filesep date filesep 'stick']);
[new_epdirs.dir] = deal(date);
epdirs = [epdirs; new_epdirs];
epdirs(arrayfun(@(e) e.name(1) == '.', epdirs)) = [];
date = '20160714'; % we re-did the rest of the 3/4 trials
new_epdirs = dir([datadir filesep date filesep 'stick']);
[new_epdirs.dir] = deal(date);
epdirs = [epdirs; new_epdirs];
epdirs(arrayfun(@(e) e.name(1) == '.', epdirs)) = [];

N = length(unique(arrayfun(@(e) e.name, epdirs, 'uniformoutput',false)));
eps = struct('endeff', cell(N,1), ...
             'material', cell(N,1), ...
             'flow', cell(N,1), ...
             'data', cell(N,1), ...
             'features', cell(N,1));

for i=1:length(epdirs)
    choose_stick_4
end

%save stickdata eps

