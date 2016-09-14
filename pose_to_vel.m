function [vel, acc, posefilt, forcefilt, anglefilt] = pose_to_vel(pose, force, endeffradius)

    pose(isnan(pose)) = 0;
    force(isnan(force)) = 0;

    % unity-gain first-order low-pass filter
    a = .02;
    b = [1 .02-1];
    
    % calculate velocity = dx/dt + cross(w, r)
    %   - dx/dt at end-effector (ball center)
    %   - cross(w, r) where w = angular velocity vector at ball center, r = ball radius
    posefilt = filtfilt(a, b, pose(:,2:4));
    forcefilt = filtfilt(a, b, force(:,2:4)); % why not
    anglefilt = filtfilt(a, b, pose(:,5:end)); %for the angles
    dx = diff(posefilt);
    dt = diff(pose(:,1));
    vel = filtfilt(a, b, bsxfun(@rdivide, dx, dt));
    if nargin == 3
        for i=1:size(vel,1)
            vel(i,:) = vel(i,:) + real(logm(xfconv(pose(i+1,5:7)) \ xfconv(pose(i,5:7)))*3000 * [0 0 -endeffradius]')';
        end
    end
    vel = [vel; vel(end,:)]; % duplicate the last point to make it the same length
    acc = filtfilt(a, b, bsxfun(@rdivide, diff(vel), dt));
    acc = [acc; acc(end,:)];
    
end
