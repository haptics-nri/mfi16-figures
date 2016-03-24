function [vel, acc, posefilt, forcefilt] = pose_to_vel(pose, force)

    % unity-gain first-order low-pass filter
    a = .02;
    b = [1 .02-1];
    
    % calculate velocity = dx/dt + cross(w, r)
    %   - dx/dt at end-effector (ball center)
    %   - cross(w, r) where w = angular velocity vector at ball center, r = ball radius
    posefilt = filtfilt(a, b, pose(:,2:4));
    forcefilt = filtfilt(a, b, force(:,2:4)); % why not
    dx = diff(posefilt);
    dt = diff(pose(:,1));
    vel = filtfilt(a, b, bsxfun(@rdivide, dx, dt));
    vel = [vel; vel(end,:)]; % duplicate the last point to make it the same length
    acc = filtfilt(a, b, bsxfun(@rdivide, diff(vel), dt));
    acc = [acc; acc(end,:)];
    
end
