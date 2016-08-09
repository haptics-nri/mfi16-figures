function sl = visforce(force, pose, video, video_offset, mass, com, trialDate, trialName,vibration)
%VISFORCE Creates the graphs/visualization of the given data
%
%force is force data, pose is position data
%video is the name/location of the video file
%video_offset is the offset between the video and the data
%mass is the mass of the end effector
%com is the center of mass of the end effector
%trialDate is the date the trial was taken
%trialName is the name of the trial
%vibration is the vibration data


addpath(genpath('RANSAC-Toolbox'))
addpath(genpath('geom3d'))
addpath('vlc-matlab')
vlc_setup
%vlc_wrapper cleanup % FIXME crashes matlab sometimes

%Aligning and filtering data
[start, stop, taps] = narrow_to_taps(force);
[vel, acc, posefilt, forcefilt, anglefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));

%Making initial position point (0,0,0)
posefilt(:,1) = posefilt(:,1) - posefilt(1,1);
posefilt(:,2) = posefilt(:,2) - posefilt(1,2);
posefilt(:,3) = posefilt(:,3) - posefilt(1,3);


% correct for non-quasistatic motion
forcefilt = forcefilt - mass*acc/1000;

%Various ranges/scaling factors
xyzrange = [min(posefilt(:,1)) max(posefilt(:,1)) ...
    min(posefilt(:,2)) max(posefilt(:,2)) ...
    min(posefilt(:,3)) max(posefilt(:,3))];
xyzscale = max(xyzrange([2 4 6]) - xyzrange([1 3 5]))/2;
xyfrange = [min(min(forcefilt(:,1:2))) max(max(forcefilt(:,1:2)))];
xyfscale = xyfrange(2) - xyfrange(1);
zfrange  = [min(forcefilt(:,3)) max(forcefilt(:,3))];
zfscale  = zfrange(2) - zfrange(1);
vrange   = [min(min(vel(:,1:2))) max(max(vel(:,1:2)))];
vscale   = vrange(2) - vrange(1);
arange   = [min(min(acc(:,1:2))) max(max(acc(:,1:2)))];
ascale   = arange(2) - arange(1);

%Equations to create lines later
scale = 1.7;
xyf_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + forcefilt(i-start+1,1)*xyzscale/xyfscale * scale];
xyf_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + forcefilt(i-start+1,2)*xyzscale/xyfscale * scale];
xyf_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
zf_x  = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1)];
zf_y  = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2)];
zf_z  = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3) + forcefilt(i-start+1,3)*xyzscale/zfscale * scale];
vel_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + vel(i-start+1,1)*xyzscale/vscale * scale];
vel_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + vel(i-start+1,2)*xyzscale/vscale * scale];
vel_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
acc_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + acc(i-start+1,1)*xyzscale/ascale * scale ];
acc_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + acc(i-start+1,2)*xyzscale/ascale * scale];
acc_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];

%cyl_xyz = @(i) nth(6, [1 2 3], @cylinder_surf, 1, 5, 0, xfconv(pose(i,5:7))*[0 0 -1]', posefilt(i-start+1,1:3), 1, [0 0 0]);

txt_x = @(i) posefilt(i-start+1,1) + 2;
txt_y = @(i) posefilt(i-start+1,2) + 2;
txt_z = @(i) posefilt(i-start+1,3) - 3;
% txt_text = @(i) sprintf('Position: (%g, %g, %g) mm\nVelocity: %g mm/s\nForce: (%g, %g) N\n\\mu=%g', ...
%     posefilt(i-start+1,1), posefilt(i-start+1,2), posefilt(i-start+1,3), ...
%     sqrt(sum(vel(i-start+1,1:2).^2,2)), ...
%     sqrt(sum(forcefilt(i-start+1,1:2).^2)),  forcefilt(i-start+1,3), ...
%     sqrt(sum(forcefilt(i-start+1,1:2).^2)) / forcefilt(i-start+1,3));

%Function to get values for table
txt_text = @(i) {sprintf('%g mm', posefilt(i-start+1,1)); ...
    sprintf('%g mm',posefilt(i-start+1,2)); ...
    sprintf('%g mm',posefilt(i-start+1,3)); ...
    sprintf('%g mm/s',sqrt(sum(vel(i-start+1,1:2).^2,2))); ...
    sprintf('%g N',sqrt(sum(forcefilt(i-start+1,1:2).^2))); ...
    sprintf('%g N',forcefilt(i-start+1,3)); ...
    sprintf('%g', sqrt(sum(forcefilt(i-start+1,1:2).^2)) / forcefilt(i-start+1,3))};


%Slider functions
pc_x = @(i) [pose(i,1) pose(i,1)] - pose(1,1);
fc_x = @(i) [force(i,1) force(i,1)] - force(1,1);
vib_x = @(i) [vibration(i,1) vibration(i,1)] - force(1,1);

%Function for title text
title_text = @(i) sprintf('t=%g s (i=%d)', force(i,1) - force(1,1), i);

clf;
%Puts name of trial on top of the figure header
set(gcf,'NumberTitle','off','Name', [trialName ': ' trialDate])

%General layout of figure
p = get(gcf, 'Position');
set(gcf, 'Position', [p(1:2) 1024 512]);
ax_slider = subplot(2,3,[1 4]);
ax_pose   = subplot(3,2,2); 
ax_force  = subplot(3,2,4);
ax_vibrate = subplot(3,2,6);
%set(gcf, 'CloseRequestFcn', @(s,c) {funwrap(@vlc_wrapper, {'cleanup'}) close});
% fucking matlab crashes all the fucking time when I try to do this ^

%Plotting/updating the plot for the slider
axes(ax_slider);
hold on;

    %Making the cylinder to show position/rotation
    [xcyl, ycyl, zcyl] = cylinder(4);
    zcyl = zcyl*40;
    oldCylValues = {xcyl; ycyl; zcyl};
    
    %Updating and plotting the cylinder in its first position
    newCylPoints = rotateBody(oldCylValues, posefilt(1,:), anglefilt(1,:),1);
    cylPlot = surf(newCylPoints{1}, newCylPoints{2}, newCylPoints{3});
    
    %Making the sphere to show position/rotation
    [xsphere, ysphere, zsphere] = sphere;
    xsphere = xsphere * 12; ysphere = ysphere * 12; zsphere = zsphere * 12;
    oldSphereValues = {xsphere; ysphere; zsphere};
    
    %Updating and plotting the sphere in its first position
    newSpherePoints = rotateBody(oldSphereValues, posefilt(1,:), anglefilt(1,:), 2);
    spherePlot = surf(newSpherePoints{1}, newSpherePoints{2}, newSpherePoints{3});
    
%Plotting the lines (force, velocity, etc.)
xyf = line(xyf_x(start), xyf_y(start), xyf_z(start), 'Color','r','Linewidth',2);
zf  = line(zf_x(start),  zf_y(start),  zf_z(start),  'Color','r','Linewidth',2);
vel = line(vel_x(start), vel_y(start), vel_z(start), 'Color','g','Linewidth',2);
acc = line(acc_x(start), acc_y(start), acc_z(start), 'Color','b','LineWidth',2);

%Making graph look nice
hold off;
grid on;
axis equal vis3d;
siz = xyzscale + 5;
axis([xyzrange(1)-siz xyzrange(2)+siz ...
    xyzrange(3)-siz xyzrange(4)+siz ...
    xyzrange(5)-siz xyzrange(6)+siz]);
xlabel X;
ylabel Y;
zlabel Z;

%Adjusting the position of the plot slightly
slider_pos = get(ax_slider,'Position');
set(ax_slider,'Position',[slider_pos(1) slider_pos(2)+.04 slider_pos(3:4)*.93])

%Plotting/updating the plot for the position
axes(ax_pose);
plot(pose(start:stop,1)-pose(1,1), posefilt); 
a = axis;
posecur = line(pc_x(start), [a(3) a(4)], 'Color','g');
ylabel('Position (mm)');
legend('X','Y', 'Z');
pose_pos = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[pose_pos(1)+.09 pose_pos(2) pose_pos(3) pose_pos(4)]) %(1) + .04 (3) * 1.1

%Plotting/updating the plot for force
axes(ax_force);
plot(force(start:stop,1)-force(1,1), forcefilt)
a = axis;
forcecur = line(fc_x(start), [a(3) a(4)], 'Color','r');
ylabel('Force (N)');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
force_pos = get(gca,'Position'); %Changing position of plot slightly
set(gca,'Position',[force_pos(1)+.09 force_pos(2)+.03 force_pos(3) force_pos(4)])

%Plotting/updating the plot for vibration
axes(ax_vibrate)
plot(vibration(start:stop,1)-vibration(1,1),vibration(start:stop,2:end))
a = axis;
vibcur = line(vib_x(start), [a(3) a(4)], 'Color','b');
ylabel('Vibrations (m/s^2)');
legend('X','Y','Z')
xlabel('Time (s)')
vib_pos = get(gca,'Position');
set(gca,'Position',[vib_pos(1) + .09 vib_pos(2) + .06 vib_pos(3) vib_pos(4)])

[~, htitle] = suplabel(title_text(start), 't');

%Making a table for the data
tableData = txt_text(start);
%{'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
tTable = uitable('Data',tableData,'ColumnName',{''},'RowName',...
    {'X Position','Y Position','Z Position','Velocity','Tangent Force','Normal Force','mu'},...
    'Tag','Table','TooltipString','Table','Parent',gcf,'Position',[420 300 209 141]);


%Showing the video (if it exists)
if exist(video,'file')
    vh = vlc_wrapper('init');
    vp = vlc_wrapper('open', vh, video);
    vlc_wrapper('frame', vp, (video_offset+0.5-(force(round(taps(1)),1)-force(1,1)))/29.97*1000);
end

%TODO: This is a very silly way of doing this... Find a more concise
%way
if exist(video,'file')
    q = @(i) {
        set(htitle, 'String', title_text(i))
        
        %Updating the cone
        feval(@(xyz) ...
       set(cylPlot, 'Xdata',xyz{1}, ...
       'YData',xyz{2}, ...
       'ZData', xyz{3}), ...
       rotateBody(oldCylValues, posefilt(i-start+1,:), anglefilt(i-start+1,:),1))
       
        %Updating the sphere
       feval(@(xyz) ...
       set(spherePlot, 'Xdata', xyz{1}, ...
       'YData', xyz{2}, ...
       'Zdata', xyz{3}), ...
       rotateBody(oldSphereValues, posefilt(i-start+1,:), anglefilt(i-start+1,:),2))
       
        set(xyf, 'XData', xyf_x(i), ...
        'YData', xyf_y(i), ...
        'ZData', xyf_z(i))
        set(zf,  'XData', zf_x(i), ...
        'YData', zf_y(i), ...
        'ZData', zf_z(i))
        set(vel, 'XData', vel_x(i), ...
        'YData', vel_y(i), ...
        'ZData', vel_z(i))
        set(acc, 'XData', acc_x(i), ...
        'YData', acc_y(i), ...
        'ZData', acc_z(i))
%         set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
%         'String', txt_text(i))
        set(tTable,'Data',txt_text(i))


        % lockstep figure
        set(posecur,  'XData', pc_x(i))
        set(forcecur, 'XData', fc_x(i))
        set(vibcur,'XData',vib_x(i))
        
        vlc_wrapper('frame', vp, 1000*(force(i,1)-force(1,1) + (video_offset+0.5)/29.97 - (force(round(taps(1)),1)-force(1,1))))
        
        };
else
    q = @(i) {
        set(htitle, 'String', title_text(i))
        
       %Updating the cone
       feval(@(xyz) ...
       set(cylPlot, 'Xdata',xyz{1}, ...
       'YData',xyz{2}, ...
       'ZData', xyz{3}), ...
       rotateBody(oldCylValues, posefilt(i-start+1,:), anglefilt(i-start+1,:),1))
       
       %Updating the sphere
       feval(@(xyz) ...
       set(spherePlot, 'Xdata', xyz{1}, ...
       'YData', xyz{2}, ...
       'Zdata', xyz{3}), ...
       rotateBody(oldSphereValues, posefilt(i-start+1,:), anglefilt(i-start+1,:),2))
       
       %Updating various lines
        set(xyf, 'XData', xyf_x(i), ...
        'YData', xyf_y(i), ...
        'ZData', xyf_z(i))
        set(zf,  'XData', zf_x(i), ...
        'YData', zf_y(i), ...
        'ZData', zf_z(i))
        set(vel, 'XData', vel_x(i), ...
        'YData', vel_y(i), ...
        'ZData', vel_z(i))
        set(acc, 'XData', acc_x(i), ...
        'YData', acc_y(i), ...
        'ZData', acc_z(i))

        %Updating the table
         set(tTable,'Data',txt_text(i))


        % lockstep figure
        set(posecur,  'XData', pc_x(i))
        set(forcecur, 'XData', fc_x(i))
        set(vibcur,'XData',vib_x(i))
        };
    
end

%Setting the function for the slider
p = @(s,e) q(round(s.Value));

%Slider
sl = uicontrol('Style','slider', ...
    'Min',start, 'Max',stop, 'Value',start, ...
    'SliderStep',[.001 .1], ...
    'Position',[100 10 400 20], ...
    'CreateFcn',p, 'Callback',p);
uicontrol('Style','text', ...
    'Position',[45 10 50 20], ...
    'String',sprintf('%g', pose(start,1)-pose(1,1)));
uicontrol('Style','text', ...
    'Position',[505 10 50 20], ...
    'String',sprintf('%g', pose(stop,1)-pose(1,1)));


%Quit button
uicontrol('Style','pushbutton',...
    'Position', [850 10 100 40], ...
    'String', 'Quit', 'Callback',@quitfunc);

%Choose a new data set
uicontrol('Style','pushbutton',...
    'Position', [720 10 100 40], ...
    'String', 'Choose New Data','Callback',@NewDataSet) %switchData)

%Looking at the material properties
uicontrol('Style','pushbutton',...
    'Position', [590 10 100 40],...
    'String','Trial Properties','Callback',@propMenu)

    function quitfunc(~,~)
        %Closes out of the movie and the graphs
        close(gcf)
    end

    function NewDataSet(~,~)
        go_visforce
    end

    function propMenu(~,~)
        %Getting information of new/old figure windows
        propFig = figure;
        pos_propFig = get(propFig, 'Position');
        
        %Changing the size of the new figure window
        set(propFig,'Position', [pos_propFig(1:2) 500 400])
        
        foldername = video(1:(end-15));
        %Uploading .flow file for material and tooling ball (if applicable)
        if exist([foldername 'stick.flow'],'file')
            %To find out what "questions" to ask: flow.answers.keys
            flow = parse_flow([foldername 'stick.flow']);
            toolingBall = flow.answers('tooling ball diameter').text;
            toolingBall = [toolingBall ''''''];
            materialName = flow.answers('surface name').text;
        else materialName = 'Unknown';
            toolingBall = 'Unknown';
        end
        
        stf = 90;
        %Placing information on the figure window
        uicontrol('Style','Text','String','Date:', 'Position',[-10 370-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',trialDate,'Position',[180 430-stf 200 50])
        uicontrol('Style','Text','String','Trial Name:','Position',[-10 300-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',trialName,'Position',[180 360-stf 200 50])
        uicontrol('Style','Text','String','Material Name:','Position',[-10 230-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String', materialName,'Position',[180 290-stf 200 50])
        uicontrol('Style','Text','String','Tooling Ball Diam.:','Position',[-10 160-stf 200 100],...
            'FontSize',20)
        uicontrol('Style','edit','String',toolingBall,'Position',[180 220-stf 200 50])
        uicontrol('Style','Text','String','File Directory:','Position',[-10 90-stf 200 100],...
            'FontSize',20)
        hEdit = uicontrol('Style','edit','String',foldername,'Position',[180 150-stf 200 50]);
        uicontrol('Style','pushbutton','String','Close','Position',[395 5 100 50],...
            'Callback',@quitfunc)
        hCMenu = uicontextmenu;                       %# Create a context menu
        uimenu(hCMenu,'Label','Copy',...              %# Create a menu item
            'Callback',@(hObject,eventData) clipboard('copy',get(hEdit,'String')));
        set(hEdit,'UIContextMenu',hCMenu);            %# Add context menu to control
    end
end