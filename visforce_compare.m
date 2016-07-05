function [sl, sl2] = visforce_compare(force1, pose1, video1, video_offset1,trialDate1, trialName1, vibration1,...
    force2, pose2, video2, video_offset2, trialDate2, trialName2, vibration2, ...
    mass, com)

%Adding various paths
addpath(genpath('RANSAC-Toolbox'))
addpath(genpath('geom3d'))
addpath('vlc-matlab')
vlc_setup
%vlc_wrapper cleanup % FIXME crashes matlab sometimes

%Running processing stuff for the first set of data
[start1, stop1, taps1] = narrow_to_taps(force1);
[vel1, acc1, posefilt1, forcefilt1] = pose_to_vel(pose1(start1:stop1,:), force1(start1:stop1,:));

%Making initial point of first set of data (0,0,0)
posefilt1(:,1) = posefilt1(:,1) - posefilt1(1,1);
posefilt1(:,2) = posefilt1(:,2) - posefilt1(1,2);
posefilt1(:,3) = posefilt1(:,3) - posefilt1(1,3);

%Running processing stuff for the second set of data
[start2, stop2, taps2] = narrow_to_taps(force2);
[vel2, acc2, posefilt2, forcefilt2] = pose_to_vel(pose2(start2:stop2,:), force2(start2:stop2,:));

%Making initial point of first set of data (0,0,0)
posefilt2(:,1) = posefilt2(:,1) - posefilt2(1,1);
posefilt2(:,2) = posefilt2(:,2) - posefilt2(1,2);
posefilt2(:,3) = posefilt2(:,3) - posefilt2(1,3);

% correct for non-quasistatic motion for first set of data
forcefilt1 = forcefilt1 - mass*acc1/1000;

xyzrange1 = [min(posefilt1(:,1)) max(posefilt1(:,1)) ...
    min(posefilt1(:,2)) max(posefilt1(:,2)) ...
    min(posefilt1(:,3)) max(posefilt1(:,3))];
xyzscale1 = max(xyzrange1([2 4 6]) - xyzrange1([1 3 5]))/2;
xyfrange1 = [min(min(forcefilt1(:,1:2))) max(max(forcefilt1(:,1:2)))];
xyfscale1 = xyfrange1(2) - xyfrange1(1);
zfrange1  = [min(forcefilt1(:,3)) max(forcefilt1(:,3))];
zfscale1  = zfrange1(2) - zfrange1(1);
vrange1   = [min(min(vel1(:,1:2))) max(max(vel1(:,1:2)))];
vscale1   = vrange1(2) - vrange1(1);
arange1   = [min(min(acc1(:,1:2))) max(max(acc1(:,1:2)))];
ascale1   = arange1(2) - arange1(1);

% correct for non-quasistatic motion for second set of data
forcefilt2 = forcefilt2 - mass*acc2/1000;

xyzrange2 = [min(posefilt2(:,1)) max(posefilt2(:,1)) ...
    min(posefilt2(:,2)) max(posefilt2(:,2)) ...
    min(posefilt2(:,3)) max(posefilt2(:,3))];
xyzscale2 = max(xyzrange2([2 4 6]) - xyzrange2([1 3 5]))/2;
xyfrange2 = [min(min(forcefilt2(:,1:2))) max(max(forcefilt2(:,1:2)))];
xyfscale2 = xyfrange2(2) - xyfrange2(1);
zfrange2  = [min(forcefilt2(:,3)) max(forcefilt2(:,3))];
zfscale2  = zfrange2(2) - zfrange2(1);
vrange2   = [min(min(vel2(:,1:2))) max(max(vel2(:,1:2)))];
vscale2   = vrange2(2) - vrange2(1);
arange2   = [min(min(acc2(:,1:2))) max(max(acc2(:,1:2)))];
ascale2   = arange2(2) - arange2(1);

%Various anon functions needed. A lot more inputs so it works for both sets
%   of data
xyf_x = @(i,posefilt,start,forcefilt, xyzscale,xyfscale) ...
    [posefilt(i-start+1,1), posefilt(i-start+1,1) + forcefilt(i-start+1,1)*xyzscale/xyfscale];
xyf_y = @(i,posefilt,start,forcefilt,xyzscale,xyfscale) ...
    [posefilt(i-start+1,2), posefilt(i-start+1,2) + forcefilt(i-start+1,2)*xyzscale/xyfscale];
xyf_z = @(i,posefilt,start) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
zf_x  = @(i,posefilt,start) [posefilt(i-start+1,1), posefilt(i-start+1,1)];
zf_y  = @(i,posefilt,start) [posefilt(i-start+1,2), posefilt(i-start+1,2)];
zf_z  = @(i,posefilt,start,forcefilt,xyzscale,zfscale)...
    [posefilt(i-start+1,3), posefilt(i-start+1,3) + forcefilt(i-start+1,3)*xyzscale/zfscale];
vel_x = @(i,posefilt,start,vel,xyzscale,vscale)...
    [posefilt(i-start+1,1), posefilt(i-start+1,1) + vel(i-start+1,1)*xyzscale/vscale];
vel_y = @(i,posefilt,start,vel,xyzscale,vscale)...
    [posefilt(i-start+1,2), posefilt(i-start+1,2) + vel(i-start+1,2)*xyzscale/vscale];
vel_z = @(i,posefilt,start) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
acc_x = @(i,posefilt,start,acc,xyzscale,ascale)...
    [posefilt(i-start+1,1), posefilt(i-start+1,1) + acc(i-start+1,1)*xyzscale/ascale];
acc_y = @(i,posefilt,start,acc,xyzscale,ascale) [posefilt(i-start+1,2), posefilt(i-start+1,2) + acc(i-start+1,2)*xyzscale/ascale];
acc_z = @(i,posefilt,start) [posefilt(i-start+1,3), posefilt(i-start+1,3)];

cyl_xyz = @(i,pose,posefilt,start) nth(6, [1 2 3], @cylinder_surf, 1, 5, 0, xfconv(pose(i,5:7))*[0 0 -1]', posefilt(i-start+1,1:3), 1, [0 0 0]);

txt_text = @(i,posefilt,start,vel,forcefilt) {sprintf('%g mm', posefilt(i-start+1,1)); ...
    sprintf('%g mm',posefilt(i-start+1,2)); ...
    sprintf('%g mm',posefilt(i-start+1,3)); ...
    sprintf('%g mm/s',sqrt(sum(vel(i-start+1,1:2).^2,2))); ...
    sprintf('%g N',sqrt(sum(forcefilt(i-start+1,1:2).^2))); ...
    sprintf('%g N',forcefilt(i-start+1,3)); ...
    sprintf('%g', sqrt(sum(forcefilt(i-start+1,1:2).^2)) / forcefilt(i-start+1,3))};

%Slider functions (once again, more inputs so it works for both sets of
%   data)
pc_x = @(i,pose) [pose(i,1) pose(i,1)] - pose(1,1);
fc_x = @(i,force) [force(i,1) force(i,1)] - force(1,1);
vib_x = @(i,vibration,force) [vibration(i,1) vibration(i,1)] - force(1,1);

title_text = @(i,force) sprintf('t=%g s (i=%d)', force(i,1) - force(1,1), i);


clf;
%Puts name of trials on top of the figure header
set(gcf,'NumberTitle','off','Name', [trialName1 ': ' trialDate1 ' & ' trialName2 ': ' trialDate2])
set(gcf, 'Position', [85 275 1280 530]);

%Axes for first trial
ax_slider_1 = subplot(4,4,[9 10 13 14]);
ax_pose_1   = subplot(4,4,1);
ax_force_1  = subplot(4,4,5);
ax_vibrate_1 = subplot(4,4,2);

%Axes for second trial
ax_slider_2 = subplot(4,4,[11 12 15 16]);
ax_pose_2   = subplot(4,4,4);
ax_force_2  = subplot(4,4,8);
ax_vibrate_2 = subplot(4,4,3);

%%%
%Plotting/updating the plot for the slider for the first trial
axes(ax_slider_1);
hold on;
cyl1 = feval(@(xyz) surf(xyz{1}, xyz{2}, xyz{3}), cyl_xyz(start1,pose1,posefilt1,start1));
xyf1 = line(xyf_x(start1,posefilt1,start1,forcefilt1, xyzscale1,xyfscale1), ...
    xyf_y(start1,posefilt1,start1,forcefilt1,xyzscale1,xyfscale1), ...
    xyf_z(start1,posefilt1,start1), 'Color','r');
zf1  = line(zf_x(start1,posefilt1,start1),  zf_y(start1,posefilt1,start1),  ...
    zf_z(start1,posefilt1,start1,forcefilt1,xyzscale1,zfscale1),  'Color','r');
velline1 = line(vel_x(start1,posefilt1,start1,vel1,xyzscale1,vscale1),...
    vel_y(start1,posefilt1,start1,vel1,xyzscale1,vscale1),...
    vel_z(start1,posefilt1,start1), 'Color','g');
accline1 = line(acc_x(start1,posefilt1,start1,acc1,xyzscale1,ascale1),...
    acc_y(start1,posefilt1,start1,acc1,xyzscale1,ascale1),...
    acc_z(start1,posefilt1,start1), 'Color','b');

hold off;
grid on;
axis vis3d;
siz1 = xyzscale1 + 5;
axis([xyzrange1(1)-siz1 xyzrange1(2)+siz1 ...
    xyzrange1(3)-siz1 xyzrange1(4)+siz1 ...
    xyzrange1(5)-siz1 xyzrange1(6)+siz1]);
xlabel X;
ylabel Y;
zlabel Z;
%Adjusting the position of the plot slightly
slider_pos_1 = get(ax_slider_1,'Position');
set(ax_slider_1,'Position',[slider_pos_1(1)-.08 slider_pos_1(2) slider_pos_1(3:4)])
%%%

%%%
%Plotting/updating the plot for the slider for the first trial
axes(ax_slider_2);
hold on;
cyl2 = feval(@(xyz) surf(xyz{1}, xyz{2}, xyz{3}), cyl_xyz(start2,pose2,posefilt2,start2));
xyf2 = line(xyf_x(start2,posefilt2,start2,forcefilt2, xyzscale2,xyfscale2), ...
    xyf_y(start2,posefilt2,start2,forcefilt2,xyzscale2,xyfscale2), ...
    xyf_z(start2,posefilt2,start2), 'Color','r');
zf2  = line(zf_x(start2,posefilt2,start2),  zf_y(start2,posefilt2,start2),  ...
    zf_z(start2,posefilt2,start2,forcefilt2,xyzscale2,zfscale2),  'Color','r');
velline2 = line(vel_x(start2,posefilt2,start2,vel2,xyzscale2,vscale2),...
    vel_y(start2,posefilt2,start2,vel2,xyzscale2,vscale2),...
    vel_z(start2,posefilt2,start2), 'Color','g');
accline2 = line(acc_x(start2,posefilt2,start2,acc2,xyzscale2,ascale2),...
    acc_y(start2,posefilt2,start2,acc2,xyzscale2,ascale2),...
    acc_z(start2,posefilt2,start2), 'Color','b');

hold off;
grid on;
axis vis3d;
siz2 = xyzscale2 + 5;
axis([xyzrange2(1)-siz2 xyzrange2(2)+siz2 ...
    xyzrange2(3)-siz2 xyzrange2(4)+siz2 ...
    xyzrange2(5)-siz2 xyzrange2(6)+siz2]);
xlabel X;
ylabel Y;
zlabel Z;
%Adjusting the position of the plot slightly
slider_pos_2 = get(ax_slider_2,'Position');
set(ax_slider_2,'Position',[slider_pos_2(1)+.08 slider_pos_2(2) slider_pos_2(3:4)])
%%%



%Plotting/updating the plot for the position for the first trial
axes(ax_pose_1);
plot(pose1(start1:stop1,1)-pose1(1,1), posefilt1);
a = axis;
posecur1 = line(pc_x(start1,pose1), [a(3) a(4)], 'Color','g');
ylabel('Position (mm)');
legend('X','Y', 'Z');
pose_pos_1 = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[pose_pos_1(1)-.08 pose_pos_1(2) pose_pos_1(3) pose_pos_1(4)]) %(1) + .04 (3) * 1.1

%Plotting/updating the plot for the position for the second trial
axes(ax_pose_2);
plot(pose2(start2:stop2,1)-pose2(1,1), posefilt2);
a = axis;
posecur2 = line(pc_x(start2,pose2), [a(3) a(4)], 'Color','g');
ylabel('Position (mm)');
legend('X','Y', 'Z');
pose_pos_2 = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[pose_pos_2(1)+.08 pose_pos_2(2) pose_pos_2(3) pose_pos_2(4)]) %(1) + .04 (3) * 1.1

%Plotting/updating the plot for force for the first trial
axes(ax_force_1);
plot(force1(start1:stop1,1)-force1(1,1), forcefilt1)
a = axis;
forcecur1 = line(fc_x(start1,force1), [a(3) a(4)], 'Color','r');
ylabel('Force (N)');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
force_pos_1 = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[force_pos_1(1)-.08 force_pos_1(2) force_pos_1(3) force_pos_1(4)])

%Plotting/updating the plot for force for the second trial
axes(ax_force_2);
plot(force2(start2:stop2,1)-force2(1,1), forcefilt2)
a = axis;
forcecur2 = line(fc_x(start2,force2), [a(3) a(4)], 'Color','r');
ylabel('Force (N)');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
force_pos_2 = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[force_pos_2(1)+.08 force_pos_2(2) force_pos_2(3) force_pos_2(4)])

%Plotting/updating the plot for vibration for the first trial
axes(ax_vibrate_1)
plot(vibration1(start1:stop1,1)-vibration1(1,1),vibration1(start1:stop1,2:end))
a = axis;
vibcur1 = line(vib_x(start1,vibration1,force1), [a(3) a(4)], 'Color','b');
ylabel('Vibrations (m/s^2)');
legend('X','Y','Z')
xlabel('Time (s)')
vib_pos_1 = get(gca,'Position');
set(gca,'Position',[vib_pos_1(1)-.08 vib_pos_1(2) vib_pos_1(3) vib_pos_1(4)])

%Plotting/updating the plot for vibration for the second trial
axes(ax_vibrate_2)
plot(vibration2(start2:stop2,1)-vibration2(1,1),vibration2(start2:stop2,2:end))
a = axis;
vibcur2 = line(vib_x(start2,vibration2,force2), [a(3) a(4)], 'Color','b');
ylabel('Vibrations (m/s^2)');
legend('X','Y','Z')
xlabel('Time (s)')
vib_pos_2 = get(gca,'Position');
set(gca,'Position',[vib_pos_2(1)+.08 vib_pos_2(2) vib_pos_2(3) vib_pos_2(4)])


%Making a table for the data for the first trial
tableData1 = txt_text(start1,posefilt1,start1,vel1,forcefilt1);
%{'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
tTable1 = uitable('Data',tableData1,'ColumnName',{trialName1},'RowName',...
    {'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
    'Tag','Table','TooltipString','Table','Parent',gcf,'Position',[420 230 209 141],...
    'Units','Normalized');

%Making a table for the data for the second trial
tableData2 = txt_text(start2,posefilt2,start2,vel2,forcefilt2);
%{'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
tTable2 = uitable('Data',tableData2,'ColumnName',{trialName2},'RowName',...
    {'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
    'Tag','Table','TooltipString','Table','Parent',gcf,'Position',[675 230 209 141],'Units',...
    'Normalized');

%Showing the video (if it exists) for the first trial
if exist(video1,'file')
    vh1 = vlc_wrapper('init');
    vp1 = vlc_wrapper('open', vh1, video1);
    vlc_wrapper('frame', vp1, (video_offset1+0.5-(force1(round(taps1(1)),1)-force1(1,1)))/29.97*1000);
end

%Showing the video (if it exists) for the second trial
if exist(video2,'file')
    vh2 = vlc_wrapper('init');
    vp2 = vlc_wrapper('open', vh2, video2);
    vlc_wrapper('frame', vp2, (video_offset2+0.5-(force2(round(taps2(1)),1)-force2(1,1)))/29.97*1000);
end


%%%%FOR THE FIRST TRIAL
%TODO: This is a very silly way of doing this... Find a more concise
%way
if exist(video1,'file')
    q1 = @(i) {
        %set(htitle, 'String', title_text(i))
        
        % slider figure
        feval(@(xyz) ...
        set(cyl1, 'XData', xyz{1}, ...
        'YData', xyz{2}, ...
        'ZData', xyz{3}), ...
        cyl_xyz(i,pose1,posefilt1,start1))
        set(xyf1, 'XData', xyf_x(i,posefilt1,start1,forcefilt1, xyzscale1,xyfscale1), ...
        'YData', xyf_y(i,posefilt1,start1,forcefilt1,xyzscale1,xyfscale1), ...
        'ZData', xyf_z(i,posefilt1,start1))
        set(zf1,  'XData', zf_x(i,posefilt1,start1), ...
        'YData', zf_y(i,posefilt1,start1), ...
        'ZData', zf_z(i,posefilt1,start1,forcefilt1,xyzscale1,zfscale1))
        set(velline1, 'XData', vel_x(i,posefilt1,start1,vel1,xyzscale1,vscale1), ...
        'YData', vel_y(i,posefilt1,start1,vel1,xyzscale1,vscale1), ...
        'ZData', vel_z(i,posefilt1,start1))
        set(accline1, 'XData', acc_x(i,posefilt1,start1,acc1,xyzscale1,ascale1), ...
        'YData', acc_y(i,posefilt1,start1,acc1,xyzscale1,ascale1), ...
        'ZData', acc_z(i,posefilt1,start1))
        %         set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
        %         'String', txt_text(i))
        set(tTable1,'Data',txt_text(i,posefilt1,start1,vel1,forcefilt1))
        
        
        % lockstep figure
        set(posecur1,  'XData', pc_x(i,pose1))
        set(forcecur1, 'XData', fc_x(i,force1))
        set(vibcur1,'XData',vib_x(i,vibration1,force1))
        
        vlc_wrapper('frame', vp1, 1000*(force1(i,1)-force1(1,1) + (video_offset1+0.5)/29.97 - (force1(round(taps1(1)),1)-force1(1,1))))
        
        };
else
    q1 = @(i) {
        %set(htitle, 'String', title_text(i))
        
        % slider figure
        feval(@(xyz) ...
        set(cyl1, 'XData', xyz{1}, ...
        'YData', xyz{2}, ...
        'ZData', xyz{3}), ...
        cyl_xyz(i,pose1,posefilt1,start1))
        set(xyf1, 'XData', xyf_x(i,posefilt1,start1,forcefilt1, xyzscale1,xyfscale1), ...
        'YData', xyf_y(i,posefilt1,start1,forcefilt1,xyzscale1,xyfscale1), ...
        'ZData', xyf_z(i,posefilt1,start1))
        set(zf1,  'XData', zf_x(i,posefilt1,start1), ...
        'YData', zf_y(i,posefilt1,start1), ...
        'ZData', zf_z(i,posefilt1,start1,forcefilt1,xyzscale1,zfscale1))
        set(velline1, 'XData', vel_x(i,posefilt1,start1,vel1,xyzscale1,vscale1), ...
        'YData', vel_y(i,posefilt1,start1,vel1,xyzscale1,vscale1), ...
        'ZData', vel_z(i,posefilt1,start1))
        set(accline1, 'XData', acc_x(i,posefilt1,start1,acc1,xyzscale1,ascale1), ...
        'YData', acc_y(i,posefilt1,start1,acc1,xyzscale1,ascale1), ...
        'ZData', acc_z(i,posefilt1,start1))
        %         set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
        %         'String', txt_text(i))
        set(tTable1,'Data',txt_text(i,posefilt1,start1,vel1,forcefilt1))
        
        
        % lockstep figure
        set(posecur1,  'XData', pc_x(i,pose1))
        set(forcecur1, 'XData', fc_x(i,force1))
        set(vibcur1,'XData',vib_x(i,vibration1,force1))
        
        };
end
%%%%


%%%%FOR THE SECOND TRIAL
%TODO: This is a very silly way of doing this... Find a more concise
%way
if exist(video2,'file')
    q2 = @(i) {
        %set(htitle, 'String', title_text(i))
        
        % slider figure
        feval(@(xyz) ...
        set(cyl2, 'XData', xyz{1}, ...
        'YData', xyz{2}, ...
        'ZData', xyz{3}), ...
        cyl_xyz(i,pose2,posefilt2,start2))
        set(xyf2, 'XData', xyf_x(i,posefilt2,start2,forcefilt2, xyzscale2,xyfscale2), ...
        'YData', xyf_y(i,posefilt2,start2,forcefilt2,xyzscale2,xyfscale2), ...
        'ZData', xyf_z(i,posefilt2,start2))
        set(zf2,  'XData', zf_x(i,posefilt2,start2), ...
        'YData', zf_y(i,posefilt2,start2), ...
        'ZData', zf_z(i,posefilt2,start2,forcefilt2,xyzscale2,zfscale2))
        set(velline2, 'XData', vel_x(i,posefilt2,start2,vel2,xyzscale2,vscale2), ...
        'YData', vel_y(i,posefilt2,start2,vel2,xyzscale2,vscale2), ...
        'ZData', vel_z(i,posefilt2,start2))
        set(accline2, 'XData', acc_x(i,posefilt2,start2,acc2,xyzscale2,ascale2), ...
        'YData', acc_y(i,posefilt2,start2,acc2,xyzscale2,ascale2), ...
        'ZData', acc_z(i,posefilt2,start2))
        %         set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
        %         'String', txt_text(i))
        set(tTable2,'Data',txt_text(i,posefilt2,start2,vel2,forcefilt2))
        
        
        % lockstep figure
        set(posecur2,  'XData', pc_x(i,pose2))
        set(forcecur2, 'XData', fc_x(i,force2))
        set(vibcur2,'XData',vib_x(i,vibration2,force2))
        
        vlc_wrapper('frame', vp2, 1000*(force2(i,1)-force2(1,1) + (video_offset2+0.5)/29.97 - (force2(round(taps2(1)),1)-force2(1,1))))
        
        };
else
    q2 = @(i) {
        %set(htitle, 'String', title_text(i))
        
        % slider figure
        feval(@(xyz) ...
        set(cyl2, 'XData', xyz{1}, ...
        'YData', xyz{2}, ...
        'ZData', xyz{3}), ...
        cyl_xyz(i,pose2,posefilt2,start2))
        set(xyf2, 'XData', xyf_x(i,posefilt2,start2,forcefilt2, xyzscale2,xyfscale2), ...
        'YData', xyf_y(i,posefilt2,start2,forcefilt2,xyzscale2,xyfscale2), ...
        'ZData', xyf_z(i,posefilt2,start2))
        set(zf2,  'XData', zf_x(i,posefilt2,start2), ...
        'YData', zf_y(i,posefilt2,start2), ...
        'ZData', zf_z(i,posefilt2,start2,forcefilt2,xyzscale2,zfscale2))
        set(velline2, 'XData', vel_x(i,posefilt2,start2,vel2,xyzscale2,vscale2), ...
        'YData', vel_y(i,posefilt2,start2,vel2,xyzscale2,vscale2), ...
        'ZData', vel_z(i,posefilt2,start2))
        set(accline2, 'XData', acc_x(i,posefilt2,start2,acc2,xyzscale2,ascale2), ...
        'YData', acc_y(i,posefilt2,start2,acc2,xyzscale2,ascale2), ...
        'ZData', acc_z(i,posefilt2,start2))
        %         set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
        %         'String', txt_text(i))
        set(tTable2,'Data',txt_text(i,posefilt2,start2,vel2,forcefilt2))
        
        
        % lockstep figure
        set(posecur2,  'XData', pc_x(i,pose2))
        set(forcecur2, 'XData', fc_x(i,force2))
        set(vibcur2,'XData',vib_x(i,vibration2,force2))
        };
end
%%%%

%Anon function for sliders
p1 = @(s,e) q1(round(s.Value));
p2 = @(s,e) q2(round(s.Value));

%Slider for trial 1
sl = uicontrol('Style','slider', ...
    'Min',start1, 'Max',stop1, 'Value',start1, ...
    'SliderStep',[.001 .1], ...
    'Position',[100 10 400 20], ...
    'CreateFcn',p1, 'Callback',p1);
uicontrol('Style','text', ...
    'Position',[45 10 50 20], ...
    'String',sprintf('%g', pose1(start1,1)-pose1(1,1)));
uicontrol('Style','text', ...
    'Position',[505 10 50 20], ...
    'String',sprintf('%g', pose1(stop1,1)-pose1(1,1)));

%Slider for trial 2
sl2 = uicontrol('Style','slider', ...
    'Min',start2, 'Max',stop2, 'Value',start2, ...
    'SliderStep',[.001 .1], ...
    'Position',[800 10 400 20], ...
    'CreateFcn',p2, 'Callback',p2);
uicontrol('Style','text', ...
    'Position',[745 10 50 20], ...
    'String',sprintf('%g', pose2(start2,1)-pose2(1,1)));
uicontrol('Style','text', ...
    'Position',[1205 10 50 20], ...
    'String',sprintf('%g', pose2(stop2,1)-pose2(1,1)));

%Trial properties for both trials
uicontrol('Style','pushbutton',...
    'Position',[595 180 120 48],...
    'String','Trial Properties','Units','Normalized',...
    'Callback',@propMenu)

%Change data
uicontrol('Style','pushbutton',...
    'Position',[595 125 120 48],...
    'String','New Comparison','Units','Normalized',...
    'Callback',@Compare)

%Text saying "Trial 1" and "Trial 2"
uicontrol('Style','text','String','Trial 1','Position',[195 495 200 25],...
    'FontSize',20,'Units','Normalized')
uicontrol('Style','text','String','Trial 2','Position',[925 495 200 25],...
    'FontSize',20,'Units','Normalized')

%Quit
uicontrol('Style','pushbutton',...
    'Position',[595 70 120 48],...
    'String','Quit','Units','Normalized','Callback',@quitfunc)

    function quitfunc(~,~)
        %Closes out of the movie and the graphs
        close(gcf)
    end

    function Compare(~,~)
        go_visforce_compare
    end

    function propMenu(~,~)
        %Getting information of new/old figure windows
        propFig = figure;
        pos_propFig = get(propFig, 'Position');
        
        %Changing the size of the new figure window
        set(propFig,'Position', [pos_propFig(1:2) 750 430])
        
        %Getting some information about the first trial
        foldername{1} = video1(1:(end-15));
        trialName{1} = trialName1;
        trialDate{1} = trialDate1;
        
        %Getting some information about the second trial
        foldername{2} = video2(1:(end-15));
        trialName{2} = trialName2;
        trialDate{2} = trialDate2;
        
        %Initializing
        toolingBall = cell(1,2);
        materialName = cell(1,2);
        
        %For both of the trials, figure out tooling ball size and material
        %   name
        for i = 1:2
            %Uploading .flow file for material and tooling ball (if applicable)
            if exist([foldername{i} 'stick.flow'],'file')
                %To find out what "questions" to ask: flow.answers.keys
                flow = parse_flow([foldername{i} 'stick.flow']);
                toolingBall{i} = flow.answers('tooling ball diameter').text;
                toolingBall{i} = [toolingBall{i} ''''''];
                materialName{i} = flow.answers('surface name').text;
            else materialName{i} = 'Unknown';
                toolingBall{i} = 'Unknown';
            end
        end
        stf = 90;
        
        %Placing information on the figure window
        %    Text and Edit text boxes
        uicontrol('Style','Text','String','Trial 1','Position',[180 320 200 100], ...
            'FontSize',20)
        uicontrol('Style','Text','String','Trial 2','Position',[420 320 200 100], ...
            'FontSize',20)
        uicontrol('Style','Text','String','Date:', 'Position',[-10 370-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',trialDate{1},'Position',[180 430-stf 200 50])
        uicontrol('Style','edit','String',trialDate{2},'Position',[420 430-stf 200 50])
        uicontrol('Style','Text','String','Trial Name:','Position',[-10 300-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',trialName{1},'Position',[180 360-stf 200 50])
        uicontrol('Style','edit','String',trialName{2},'Position',[420 360-stf 200 50])
        uicontrol('Style','Text','String','Material Name:','Position',[-10 230-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String', materialName{1},'Position',[180 290-stf 200 50])
        uicontrol('Style','edit','String', materialName{2},'Position',[420 290-stf 200 50])
        uicontrol('Style','Text','String','Tooling Ball Diam.:','Position',[-10 160-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',toolingBall{1},'Position',[180 220-stf 200 50])
        uicontrol('Style','edit','String',toolingBall{2},'Position',[420 220-stf 200 50])
        uicontrol('Style','Text','String','File Directory:','Position',[-10 90-stf 200 100],...
            'FontSize',20)
        hEdit = uicontrol('Style','edit','String',foldername{1},'Position',[180 150-stf 200 50]);
        hEdit2 = uicontrol('Style','edit','String',foldername{2},'Position',[420 150-stf 200 50]);
        
        %Allowing copy/pasting of the folder directories
        hCMenu = uicontextmenu;                       %# Create a context menu
        uimenu(hCMenu,'Label','Copy',...              %# Create a menu item
            'Callback',@(hObject,eventData) clipboard('copy',get(hEdit,'String')));
        set(hEdit,'UIContextMenu',hCMenu);            %# Add context menu to control
        
        hCMenu2 = uicontextmenu;
        uimenu(hCMenu2,'Label','Copy',...
            'Callback',@(hObject,eventData) clipboard('copy',get(hEdit2,'String')));
        set(hEdit2,'UIContextMenu',hCMenu2);
        
        %Quit button (closes the figure window)
        uicontrol('Style','pushbutton','String','Close','Position',[645 5 100 50],...
            'Callback',@quitfunc);
        
    end


end
