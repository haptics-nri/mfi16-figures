% part 5 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% sync vicon/force data

for i=1:length(eps)
    choose_stick_6
end
% loop through and fixup any that are too far from the mean
for idx = find(abs([eps.offset] - mean([eps.offset])) > 0.25)
    choose_stick_7
end

%save stickdata eps

