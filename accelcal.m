%% setup

addpath(genpath('RANSAC-Toolbox'));
DATADIR = '/Volumes/shared/Projects/Proton Pack/Data/20160830';
episode = 'accelcal';
number = '1';

%% load data

[~,~,da,~,~,a] = load_stick([DATADIR filesep episode filesep number filesep]);

%% fit spheres

[c1, r1] = sphereFit_ransac(a(:,2:4));
[c2, r2] = sphereFit_ransac(a(:,5:7));
[cd, rd] = sphereFit_ransac(da(:,2:4), 1e6);

%% figure

figure;
plot(a(:,1)-a(1,1), ([-1 0 0; 0 1 0; 0 0 1]  * bsxfun(@minus, a(:,2:4), c1)'/r1)', ...
     a(:,1)-a(1,1), ([-1 0 0; 0 1 0; 0 0 -1] * bsxfun(@minus, a(:,5:7), c2)'/r2)', ...
     da(:,1)-a(1,1),                           bsxfun(@minus, da(:,2:4), cd)/rd);
legend('a1-x', 'a1-y', 'a1-z', ...
       'a2-x', 'a2-y', 'a3-z', ...
       'da-x', 'da-y', 'da-z');
   
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Calibrated digital + analog acceleration');
