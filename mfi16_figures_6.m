% part 6 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% process calibrations: use d+R to calculate x+y+z for H_vic2bod + H_bal2imu

H_vic2bod = [ 0.9912 -0.0238 -0.1302    0
              0.0164  0.9982 -0.0575   36
              0.1314  0.0549  0.9898 -511.32
              0       0       0         1   ];
H_m402bod = [ 0       0       1       108.99
              1       0       0         0.53
              0       1       0        -2.98
              0       0       0         1   ];
H_bal2imu = [ 1       0       0       254.84
              0       1       0         0
              0       0       1         0
              0       0       0         1   ];

syms x y z real;
R = freecalib.R;
H1 = [R [0 y z]';      0 0 0 1];
H2 = [eye(3) [x 0 0]'; 0 0 0 1];
S = solve(H1 * H2 * [0 0 0 1]' == [d 1]');

H_vic2bod(1:3,1:3) = R;
H_vic2bod(2,4) = S.y;
H_vic2bod(3,4) = S.z;
H_bal2imu(1,4) = S.x;

