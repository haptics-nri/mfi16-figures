% given force and position/orientation, find coefficient of dynamic friction
% force is Nx4 (time, Fx, Fy, Fz)
% pose is Nx6  (time, X, Y, Z, R1, R2, R3)
% assumption: time columns are the same
% assumption: all motion is in the horizontal plane
% assumption: constant contact from last-early-spike+0.5s to first-late-spike-0.5s
function [mu_k, err] = extract_friction(force, pose)

    [start, stop] = narrow_to_taps(force);
    
    vel = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    
    % add contribution from angular velocity
    avel = vel;
    dt = diff(pose(start:stop,1));
    for i=1:size(vel,1)-1
        Rcurr = xfconv(pose(i,5:7));
        Rnext = xfconv(pose(i+1,5:7));
        avel(i,:) = vel(i,:)' + (Rnext - Rcurr)/dt(i) * Rnext' * [0 0 -3/8*25.4]';
    end
    
    unitvel = bsxfun(@rdivide, avel(:,1:2), sqrt(sum(avel(:,1:2).^2,2))); % project velocity into plane and normalize
    proj = zeros(stop-start+1,1);
    for i=1:length(proj)
        proj(i,:) = dot(-unitvel(i,:), force(start+i-1,2:3));
    end
    
    [mu_k, stats] = robustfit(force(start:stop,4), proj, [], [], 'off');
    err = stats.se;
    % TODO robust or RANSAC?
    % TODO which error metric to use
    
    % sources of error
    %   - orientation of body frame
    %        - recover from gravity calibration data
    %   - angular velocity of contact point
    % other features
    %   - histogram of contact point speed vs normal force vs vibration power (DFT321)

    %keyboard
    
end
