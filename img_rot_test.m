%%
addpath(genpath('RANSAC-Toolbox'));

%%

%[data, ~, calib] = icra17_load('../../nri/data', '20170221', 'stickcam', @(x) false, {'20170220', '2'});
%[data, ~, calib] = icra17_load('../../nri/data', '20170225', 'stickcam', @(x) false, {'20170220', '2'});
%[data, ~, calib] = icra17_load('../../nri/data', '20170304', 'stickcam', @(x) false, {'20170220', '2'});
[data, ~, calib] = icra17_load('../../nri/data', '20170303', 'stickcam', @(x) false, {'20170220', '2'});
load icra17_final H_vic2bod H_m402bod H_bal2imu
data = icra17_process('bluefox', data, calib.mass, H_vic2bod, H_m402bod, H_bal2imu);

%%
%prefix = '../../nri/data/20170221/stickcam/45';
%d = data('Wood 5');
%a = 9098;
%b = 80120;
%prefix = '../../nri/data/20170225/stickcam/4';
%d = data('Smooth wooden table');
%a = 16320;
%b = 104600;
%prefix = '../../nri/data/20170304/stickcam/5';
%d = data('Yellow felt');
%a = 53400;
%b = 23030;
prefix = '../../nri/data/20170303/stickcam/15';
d = data('Brick');
a = 19220;
b = 53900;

bt = readtable([prefix '/bluefox/bluefox_times.csv']);
[~,fa] = min(abs(d.bvei(a,1) - bt.UnixTimestamp));
[~,fb] = min(abs(d.bvei(b,1) - bt.UnixTimestamp));
fa = bt.FrameNumber(fa);
fb = bt.FrameNumber(fb);

imga = imread(sprintf([prefix '/bluefox/bluefox%d.png'], fa));
imgb = imread(sprintf([prefix '/bluefox/bluefox%d.png'], fb));

rect = [655 1200 928 704];
pt = [778 760];
imgb_mask = uint8(zeros(size(imgb)));
imgb_mask(rect(4):rect(2), rect(1):rect(3), :) = imgb(rect(4):rect(2), rect(1):rect(3), :);
imga(rect(4):rect(2), rect(1):rect(3), :) = 0;
imgb(rect(4):rect(2), rect(1):rect(3), :) = 0;

%%

apr = d.april;
apra = apr(fa);
aprb = apr(fb);

[~,ia,ib] = intersect(apra.ids, aprb.ids);
ctra = apra.centers(:,ia)';
ctrb = aprb.centers(:,ib)';
H = fitgeotrans(ctrb, ctra, 'projective');
Hinv = fitgeotrans(ctra, ctrb, 'projective');

imgb_warp = imwarp(imgb, H, 'OutputView',imref2d(size(imgb)));

imblend = imga/2 + imgb_warp/2;% + imwarp(imgb_mask, H, 'OutputView',imref2d(size(imgb)));

clf;
subplot(211);
imshow(imblend);
%for j=1:length(apra.ids)
%    text(apra.centers(1,j), apra.centers(2,j), sprintf('%d', apra.ids(j)), 'color','r');
%end

%%

pts = [];

apra = apr(fa);
for i=sort(cell2mat(apr.keys))
    aprb = apr(i);
    [~,ia,ib] = intersect(apra.ids, aprb.ids);
    if length(ia) >= 4
        ctra = apra.centers(:,ia)';
        ctrb = aprb.centers(:,ib)';
        H = fitgeotrans(ctrb, ctra, 'projective');

        ptb = [pt 1] * H.T;
        if ptb > 0
            pts(end+1,:) = ptb(1:2);
        end
    end
end

bbox = round([min(pts(:,1))-25 max(pts(:,2))+25 max(pts(:,1))+25 min(pts(:,2))-25]);

hold on;
    plot(pts(:,1), pts(:,2), 'r.', 'markersize',10);
hold off;
subplot(212);
imshow(imga(bbox(4):bbox(2), bbox(1):bbox(3), :)*3);
