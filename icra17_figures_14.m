% part 14 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% accelerometer comparison figure

[v,f,da,dg,~,a,~,dt,~,~,m] = load_stick('../../nri/data/20160906/stickcam/6/');
t = cumsum(bitand(dt, 65535))/1e6;
figure
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(211)
plot(t, bsxfun(@minus, a(:,2:4), mean(a(:,2:4))))
subplot(212)
plot(t, bsxfun(@minus, a(:,5:7), mean(a(:,5:7))))
subplot(211)
axis([9 17 -8 8])
subplot(212)
axis([9 17 -8 8])
xlabel('Time (s)')
ylabel('Internal accelerometer signal (m/s^2)')
subplot(211)
xlabel('Time (s)')
ylabel('External accelerometer signal (m/s^2)')
print -dpdf -r0 icra17_accel_compare.pdf

