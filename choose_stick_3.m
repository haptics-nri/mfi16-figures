% part 3 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% load texture data

date = '20160620';
epdirs = dir([datadir filesep date filesep 'stick']);
epdirs(arrayfun(@(e) e.name(1) == '.', epdirs)) = [];

eps = struct('endeff', cell(length(epdirs),1), ...
             'material', cell(length(epdirs),1), ...
             'flow', cell(length(epdirs),1), ...
             'data', cell(length(epdirs),1), ...
             'features', cell(length(epdirs),1));
for i=1:length(epdirs)
    choose_stick_4
end

save stickdata eps

