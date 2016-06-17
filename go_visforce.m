% This file demonstrates how to load a dataset and show the position/force
% visualizer. The dataset is from March 2016.

mass = 0.1503;
com = [-0.0029   -0.0030    0.0348 ]';
H_vic2bod = [ 0.9912   -0.0236   -0.1303         0
              0.0162    0.9982   -0.0571   36.0685
              0.1314    0.0545    0.9898 -511.6330
                   0         0         0    1.0000 ];
H_m402bod = [      0         0    1.0000  108.9900
              1.0000         0         0    0.5300
                   0    1.0000         0   -2.9800
                   0         0         0    1.0000 ];
H_bal2imu = [ 1.0000         0         0  254.3402
                   0    1.0000         0         0
                   0         0    1.0000         0
                   0         0         0    1.0000 ];
[v, f, ~,~,~, a] = load_stick('/Users/alex/Documents/research/proton/code/nri/data/20160310/black1stick/');
[~,~,~,~,~,~, v, ~,~,~,~, f] = process_stick(v, f, a, mass, com, H_vic2bod, H_m402bod, H_bal2imu, -2.5887);
visforce(f, v, '/Users/alex/Documents/research/proton/code/nri/data/20160310/black1stick/video_small.mov', 20, mass, com);
