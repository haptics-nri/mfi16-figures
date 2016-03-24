% given force and position/orientation, find "springiness" of surface
% force is Nx4 (time, Fx, Fy, Fz)
% pose is Nx6  (time, X, Y, Z, R1, R2, R3)
% vibe is Nx4  (time, Ax, Ay, Az)
% same assumptions as extract_friction and extract_springiness
function histo = extract_vibration(force, pose, vibe, mass)

    [start, stop] = narrow_to_taps(force);
    
    [vel, accfilt, ~, forcefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    speed = sqrt(sum(vel.^2,2));
    forcefiltsub = forcefilt - mass*accfilt/1000;
    forcemag = sqrt(sum(forcefiltsub.^2,2));
    
    % add contribution from angular velocity
    avel = vel;
    dt = diff(pose(start:stop,1));
    for i=1:size(vel,1)-1
        Rcurr = xfconv(pose(i,5:7));
        Rnext = xfconv(pose(i+1,5:7));
        avel(i,:) = vel(i,:)' + (Rnext - Rcurr)/dt(i) * Rnext' * [0 0 -3/8*25.4]';
    end
    
    dft = dft321(vibe(start:stop,2:4));
    
    % histogram of speed + normal force + vibration power
    nbins = 20;
    speed_bins = linspace(    0,  250, nbins+1); % mm/s
    force_bins = linspace(    0,   30, nbins+1); % N
    power_bins = linspace(-4000, 4000, nbins+1); % units?
    
    assert(~any(speed < speed_bins(1)));
    assert(~any(speed >= speed_bins(end)));
    assert(~any(forcemag < force_bins(1)));
    assert(~any(forcemag >= force_bins(end)));
    assert(~any(dft < power_bins(1)));
    assert(~any(dft >= power_bins(end)));
    
    histo = zeros(nbins, nbins, nbins);
    for s=1:nbins
        for f=1:nbins
            for p=1:nbins
                histo(s,f,p) = nnz(speed    >= speed_bins(s) & speed    < speed_bins(s+1) & ...
                                   forcemag >= force_bins(f) & forcemag < force_bins(f+1) & ...
                                   dft      >= power_bins(p) & dft      < power_bins(p+1));
            end
        end
    end
    
    assert(sum(histo(:)) == stop - start + 1);
    histo = histo / (stop - start + 1);
    
    % take into account vibration power + frequency, not just magnitude!
    %      window for frequency?
    
end
