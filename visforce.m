function sl = visforce(force, pose, video, video_offset, mass, com)

    addpath(genpath('RANSAC-Toolbox'))
    addpath(genpath('geom3d'))
    addpath('vlc-matlab')
    vlc_setup
    %vlc_wrapper cleanup % FIXME crashes matlab sometimes

    [start, stop, taps] = narrow_to_taps(force);
    [vel, acc, posefilt, forcefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    
    % correct for non-quasistatic motion
    forcefilt = forcefilt - mass*acc/1000;
    
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
    
    xyf_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + forcefilt(i-start+1,1)*xyzscale/xyfscale];
    xyf_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + forcefilt(i-start+1,2)*xyzscale/xyfscale];
    xyf_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
    zf_x  = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1)];
    zf_y  = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2)];
    zf_z  = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3) + forcefilt(i-start+1,3)*xyzscale/zfscale];
    vel_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + vel(i-start+1,1)*xyzscale/vscale];
    vel_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + vel(i-start+1,2)*xyzscale/vscale];
    vel_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
    acc_x = @(i) [posefilt(i-start+1,1), posefilt(i-start+1,1) + acc(i-start+1,1)*xyzscale/ascale];
    acc_y = @(i) [posefilt(i-start+1,2), posefilt(i-start+1,2) + acc(i-start+1,2)*xyzscale/ascale];
    acc_z = @(i) [posefilt(i-start+1,3), posefilt(i-start+1,3)];
    
    cyl_xyz = @(i) nth(6, [1 2 3], @cylinder_surf, 1, 5, 0, xfconv(pose(i,5:7))*[0 0 -1]', posefilt(i-start+1,1:3), 1, [0 0 0]);
    
    txt_x = @(i) posefilt(i-start+1,1) + 2;
    txt_y = @(i) posefilt(i-start+1,2) + 2;
    txt_z = @(i) posefilt(i-start+1,3) - 3;
    txt_text = @(i) sprintf('Position: (%g, %g, %g) mm\nVelocity: %g mm/s\nForce: (%g, %g) N\n\\mu=%g', ...
                            posefilt(i-start+1,1), posefilt(i-start+1,2), posefilt(i-start+1,3), ...
                            sqrt(sum(vel(i-start+1,1:2).^2,2)), ...
                            sqrt(sum(forcefilt(i-start+1,1:2).^2)),  forcefilt(i-start+1,3), ...
                            sqrt(sum(forcefilt(i-start+1,1:2).^2)) / forcefilt(i-start+1,3));
    
    pc_x = @(i) [pose(i,1) pose(i,1)] - pose(1,1);
    fc_x = @(i) [force(i,1) force(i,1)] - force(1,1);
    
    title_text = @(i) sprintf('t=%g s (i=%d)', force(i,1) - force(1,1), i);
    
    clf;
    p = get(gcf, 'Position');
    set(gcf, 'Position', [p(1:2) 1024 512]);
    ax_slider = subplot(2,2,[1 3]);
    ax_pose   = subplot(2,2,2);
    ax_force  = subplot(2,2,4);
    %set(gcf, 'CloseRequestFcn', @(s,c) {funwrap(@vlc_wrapper, {'cleanup'}) close});
        % fucking matlab crashes all the fucking time when I try to do this ^
    
    axes(ax_slider);
    hold on;
    cyl = feval(@(xyz) surf(xyz{1}, xyz{2}, xyz{3}), cyl_xyz(start));
    xyf = line(xyf_x(start), xyf_y(start), xyf_z(start), 'Color','r');
    zf  = line(zf_x(start),  zf_y(start),  zf_z(start),  'Color','r');
    vel = line(vel_x(start), vel_y(start), vel_z(start), 'Color','g');
    acc = line(acc_x(start), acc_y(start), acc_z(start), 'Color','b');
    txt = text(txt_x(start), txt_y(start), txt_z(start), txt_text(start));

    hold off;
    grid on;
    axis vis3d;
    siz = xyzscale + 5;
    axis([xyzrange(1)-siz xyzrange(2)+siz ...
          xyzrange(3)-siz xyzrange(4)+siz ...
          xyzrange(5)-siz xyzrange(6)+siz]);
    xlabel X;
    ylabel Y;
    zlabel Z;
    
    axes(ax_pose);
    plot(pose(start:stop,1)-pose(1,1), posefilt);
    a = axis;
    posecur = line(pc_x(start), [a(3) a(4)], 'Color','g');
    ylabel('Position (mm)');
    legend X Y Z;
    axes(ax_force);
    plot(force(start:stop,1)-force(1,1), forcefilt);
    a = axis;
    forcecur = line(fc_x(start), [a(3) a(4)], 'Color','r');
    ylabel('Force (N)');
    legend X Y Z;
    xlabel('Time (s)');
    
    [~, htitle] = suplabel(title_text(start), 't');

    vh = vlc_wrapper('init');
    vp = vlc_wrapper('open', vh, video);
    vlc_wrapper('frame', vp, (video_offset+0.5-(force(round(taps(1)),1)-force(1,1)))/29.97*1000);

    
    q = @(i) {
                 set(htitle, 'String', title_text(i))
    
                 % slider figure
                 feval(@(xyz) ...
                           set(cyl, 'XData', xyz{1}, ...
                                    'YData', xyz{2}, ...
                                    'ZData', xyz{3}), ...
                       cyl_xyz(i))
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
                 set(txt, 'Position', [txt_x(i), txt_y(i), txt_z(i)], ...
                          'String', txt_text(i))
                 
                 % lockstep figure
                 set(posecur,  'XData', pc_x(i))
                 set(forcecur, 'XData', fc_x(i))
                 
                 vlc_wrapper('frame', vp, 1000*(force(i,1)-force(1,1) + (video_offset+0.5)/29.97 - (force(round(taps(1)),1)-force(1,1))))
             };
    p = @(s,e) q(round(s.Value));
    
    sl = uicontrol('Style','slider', ...
              'Min',start, 'Max',stop, 'Value',start, ...
              'SliderStep',[.001 .1], ...
              'Position',[100 10 400 20], ...
              'CreateFcn',p, 'Callback',p);
    uicontrol('Style','text', ...
              'Position',[45 10 50 20], ...
              'String',sprintf('%g',posefilt(start,1)-posefilt(1,1)));
    uicontrol('Style','text', ...
              'Position',[505 10 50 20], ...
              'String',sprintf('%g',posefilt(end,1)-posefilt(1,1)));
end
