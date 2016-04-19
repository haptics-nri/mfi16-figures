% part 46 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% figures

clear set; % there is a var called "set" but I need to use the set() function

% sphere calibration
figure;
sphereplot(spherecalib.c, spherecalib.r, {spherecalib.x(:,2:4)});
xlabel('X position (mm)')
ylabel('Y position (mm)')
zlabel('Z position (mm)')
print -dpdf mfi16_sphere_calib.pdf

% free space calibration
figure;
plot(freecalib.int_s(:,1)-freecalib.int_s(1,1), freecalib.int_s(:,2:4))
xlabel('Time (s)')
ylabel('Measured force (N)')
legend('X component', 'Y component', 'Z component')
print -dpdf mfi16_freespace_data.pdf
figure;
transformed_grav = freecalib.ideal_grav;
for i=1:size(freecalib.ideal_grav,1)
    mfi16_figures_47
end
hold on;
h1 = plot3(freecalib.grav(:,2), freecalib.grav(:,3), freecalib.grav(:,4), '.');
h2 = plot3(transformed_grav(:,2), transformed_grav(:,3), transformed_grav(:,4), '.');
grid on
axis equal vis3d
view(11.5, 42)
xlabel X
ylabel Y
zlabel Z
legend('Measured gravity in body frame', 'World-frame gravity transformed to body frame', 'location','southeast')
print -dpdf mfi16_freespace_grav.pdf

% typical dataset
figure;
subplot(2,1,1);
plot(vbodyint{3,2,1}(:,1)-vbodyint{3,2,1}(1,1), vbodyint{3,2,1}(:,2:4));
ylabel('Position (mm)')
legend X Y Z
subplot(2,1,2);
plot(intworldsub{3,2,1}(:,1)-intworldsub{3,2,1}(1,1), intworldsub{3,2,1}(:,2:4));
xlabel('Time (s)')
ylabel('Force (N)')
legend X Y Z
print -dpdf mfi16_typical_data.pdf

