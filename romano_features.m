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
%       stmode = false (just mean) or true (mean + median + std)
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

function cells = romano_features_pre(force, pose, vibe, mass, dur, thresh, startstop)

    % 1. preprocessing

    if nargin == 7
        % taps passed in
        start = startstop(1);
        stop = startstop(2);
    else
        % find taps myself
        [start, stop] = narrow_to_taps(force);
    end
    [vel, accfilt, ~, forcefilt] = pose_to_vel(pose(start:stop,:), force(start:stop,:));
    speed = sqrt(sum(vel.^2,2));
    forcefiltsub = forcefilt - mass*accfilt/1000;
    vibe = vibe(start:stop,2:end);
    
    % 2. narrowing down to chunks
    
    chunks = zeros(0, 2);
    idx = start:stop; % between taps
    idx((speed < thresh(1)) | (abs(forcefiltsub(:,3)) < thresh(2))) = []; % speed/force thresholds
    divs = [1 find(diff(idx) > 100) length(idx)]; % divide into continuous-ish segments
    for d=1:(length(divs)-1)
        t = idx(divs(d));
        while t+floor(dur) <= idx(divs(d+1))
            chunks = [chunks
                      t-start+1 t+floor(dur)-start+1];
            t = t+floor(dur)+1;
        end
    end
    
    cells = cell(size(chunks,1), 3);
    for c=1:size(chunks,1)
        cells{c,1} = vibe(chunks(c,1):chunks(c,2), :);
        cells{c,2} = speed(chunks(c,1):chunks(c,2));
        cells{c,3} = forcefiltsub(chunks(c,1):chunks(c,2), :);
    end
    
end
   
function vectors = romano_features_post(cells, nbins, binmode, alpha, stmode)
    %persistent f;

    % 3. process each chunk
    
    %nbins = 20;
    vectors = zeros(size(cells,1), nbins + 3 + (3*stmode));
    
    for c=1:size(cells,1)
        V = cells{c,1};
        S = cells{c,2};
        F = cells{c,3};
        
        % condense and bin the vibration
        [~, V_freq] = dft321(V);
        
        switch binmode
            case 'naive'
                V_hist = zeros(1,nbins);
                N = round(length(V_freq)/nbins);
                for i=1:nbins
                    start = (i-1)*N+1;
                    stop = i*N;
                    if stop > length(V_freq)
                        if i == nbins
                            stop = length(V_freq);
                            V_hist(i) = sum(V_freq(start:stop)/(stop - start + 1));
                        else
                            V_hist(i) = 0;
                        end
                    else
                        V_hist(i) = sum(V_freq(start:stop)/(stop - start + 1));
                    end
                end
            case 'perceptual'
                V_hist = zeros(1,nbins);
                %alpha = 25;
                %if isempty(f)
                    %fprintf('regenerating f (length(V_freq) = %d)\n', length(V_freq));
                    f = linspace(0, 1500, length(V_freq));
                %end
                for i=1:nbins
                    b = i*1500/nbins;
                    w = exp(-(f - b).^2 / (2*(alpha*b)^2));
                    w = w/sum(w);
                    V_hist(i) = V_freq' * w';
                    %keyboard
                end
            otherwise
                error('Unsupported bin mode')
        end
        
        % output to feature vector
        vectors(c,:) = [V_hist stats(F(:,3), stmode) stats(S, stmode) stats(sqrt(sum(F(:,1:2).^2,2)), stmode)];
    end

end

function st = stats(arr, mode)
    if mode
        st = [mean(arr) std(arr)];
    else
        st = mean(arr);
    end
end
