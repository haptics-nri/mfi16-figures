% given force and position/orientation, find "springiness" of surface
% force is Nx4 (time, Fx, Fy, Fz)
% pose is Nx6  (time, X, Y, Z, R1, R2, R3)
% assumption: time columns are the same
% assumption: all motion is in the horizontal plane
% assumption: constant contact from last-early-spike+0.5s to first-late-spike-0.5s
function [k, z0, err] = extract_springiness(force, pose, mass)

    [start, stop] = narrow_to_taps(force);
    
    [~, accfilt, posefilt, forcefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    forcefiltsub = forcefilt - mass*accfilt/1000;
    
    [p, S] = polyfit(posefilt(:,3), forcefiltsub(:,3), 1);
    k = p(1);
    z0 = -p(2)/p(1);
    err = S.normr;
    
end
