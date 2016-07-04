% part 8 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
% process data

H_vic2bod = ...
   [0.9912   -0.0236   -0.1303         0
    0.0162    0.9982   -0.0571   36.0685
    0.1314    0.0545    0.9898 -511.6330
         0         0         0    1.0000];
    
 H_bal2imu = ...
   [1.0000         0         0  254.3402
         0    1.0000         0         0
         0         0    1.0000         0
         0         0         0    1.0000];
     
 H_m402bod = ...
   [     0         0    1.0000  108.9900
    1.0000         0         0    0.5300
         0    1.0000         0   -2.9800
         0         0         0    1.0000];

for i=1:length(eps)
    choose_stick_9
end

save stickdata eps

