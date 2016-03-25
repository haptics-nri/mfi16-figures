% produces Nx33ish feature vector loosely-following Romano & KJK 2014
% N is the number of chunks found ("dur" ms each, at least thresh(1) speed, thresh(2) normal force)
%
% works in two phases selected by 'mode'
% this is to allow tuning the hyperparameters after performing train/test split
%   pre: produces a cell array with unbinned vibration and unaveraged speed/force
%       force  = 3D force with timestamps (synchronized) (taps included)
%       pose   = 3D position/orientation with timestamps (synchronized) (taps included)
%       vibe   = 3D acceleration with timestamps (synchronized) (taps included)
%       mass   = mass of end-effector
%       dur    = duration of chunks (seconds)
%       thresh = [speed_thresh force_thresh] - no chunks will come from
%                   periods of time with speed less than speed_thresh or total force
%                   magnitude less than force_thresh
%   post: bins the vibration and averages the speed/force to produce the actual feature vector
%       cells   = the output of pre
%       nbins   = number of bins for acceleration
%       binmode = 'naive' or 'perceptual' (see paper)
%       alpha   = perceptual bin tuning thing (see paper)

function varargout = romano_features(mode, varargin)
    switch mode
        case 'pre'
            varargout{:} = romano_features_pre(varargin{:});
        case 'post'
            varargout{:} = romano_features_post(varargin{:});
        otherwise
            error('Unsupported mode')
    end
end

function cells = romano_features_pre(force, pose, vibe, mass, dur, thresh)

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
    
    cells = cell(size(chunks,1), 3);
    for c=1:size(chunks,1)
        cells{c,1} = vibe(chunks(c,1):chunks(c,2), :);
        cells{c,2} = speed(chunks(c,1):chunks(c,2));
        cells{c,3} = forcefiltsub(chunks(c,1):chunks(c,2), :);
    end
    
end
   
function vectors = romano_features_post(cells, nbins, binmode, alpha)
    % 3. process each chunk
    
    %nbins = 20;
    vectors = zeros(size(cells,1), nbins+3);
    
    for c=1:size(cells,1)
        V = cells{c,1};
        S = cells{c,2};
        F = cells{c,3};
        
        % condense and bin the vibration
        [~, V_freq] = dft321(V);
        
        switch binmode
            case 'naive'
                highfreq = 1000;
                centers = linspace(highfreq/Nbins, highfreq, Nbins+1);
                V_hist = histcounts(V_freq, centers+(mean(diff(centers))/2));
            case 'perceptual'
                V_hist = zeros(1,nbins);
                %alpha = 25;
                for i=1:nbins
                    b = i*1500/nbins;
                    f = linspace(0, 1500, length(V_freq));
                    w = (f - b).^2 / (2*alpha*b)^2;
                    V_hist(i) = dot(V_freq, w);
                end
            otherwise
                error('Unsupported bin mode')
        end
        
        % output to feature vector
        % FIXME include median+std
        vectors(c,:) = [V_hist mean(F(:,3)) mean(S) mean(sqrt(sum(F(:,1:2).^2,2)))];
    end

end
