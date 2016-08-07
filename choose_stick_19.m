% part 19 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% test/train split

endeffs = struct('features', cell(4,1), ...
                 'split', cell(4,1), ...
                 'gs', cell(4,1), ...
                 'model', cell(4,1), ...
                 'predictions', cell(4,1), ...
                 'confusion', cell(4,1), ...
                 'accuracy', cell(4,1));

for i=1:length(eps)
    choose_stick_20
end
   
for i=1:4
    choose_stick_22
end

%save endeffdata endeffs

