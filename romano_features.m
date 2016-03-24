% produces Nx33 feature vector following Romano & KJK 2014
% N is the number of chunks found ("dur" ms each, at least thresh(1) speed, thresh(2) normal force)
function vectors = romano_features(force, pose, vibe, mass, dur, thresh)

    % 1. preprocessing

    [start, stop] = narrow_to_taps(force);
    [vel, accfilt, ~, forcefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    speed = sqrt(sum(vel.^2,2));
    forcefiltsub = forcefilt - mass*accfilt/1000;
    forcemag = sqrt(sum(forcefiltsub.^2,2));
    vibe = vibe(start:stop,2:end);
    
    % 2. narrowing down to chunks
    
    chunks = zeros(0, 2);
    idx = start:stop; % between taps
    idx((speed < thresh(1)) | (forcemag < thresh(2))) = []; % speed/force thresholds
    divs = [1 find(diff(idx) > 100) length(idx)]; % divide into continuous-ish segments
    dt = mean(diff(force(:,1)));
    for d=1:(length(divs)-1)
        t = idx(divs(d));
        while t+round(dur/dt) <= idx(divs(d+1))
            chunks = [chunks
                      t-start t+round(dur/dt)-start];
            t = t+round(dur/dt)+1;
        end
    end
    
    % 3. process each chunk
    
    vectors = zeros(size(chunks,1), 33);
    
    for c=1:size(chunks,1)
        F = forcefiltsub(chunks(c,1):chunks(c,2), :);
        S = speed(chunks(c,1):chunks(c,2));
        V = vibe(chunks(c,1):chunks(c,2), :);
        
        % condense and bin the vibration
        [~, V_freq] = dft321(V);
        
        if false
            % naive binning
            highfreq = 1000;
            centers = linspace(highfreq/30, highfreq, 31);
            V_hist = histcounts(V_freq, centers+(mean(diff(centers))/2));
        else
            % perceptual binning
            V_hist = zeros(1,30);
            alpha = 25;
            for i=1:30
                b = i*1500/30;
                f = linspace(0, 1500, length(V_freq));
                w = (f - b).^2 / (2*alpha*b)^2;
                V_hist(i) = dot(V_freq, w);
            end
        end
        
        % output to feature vector
        vectors(c,:) = [V_hist mean(F(:,3)) mean(S) mean(sqrt(sum(F(:,1:2).^2,2)))];
    end

end
