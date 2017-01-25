% calculate Steinbach scan-free features
% inputs:
%   - feats: (cell array of strings) features to calculate
%       grep for strcmp to see all possibilities
%       accel features are named like the WHC paper, and sound features like the ToH paper
%   - hac, imp, mov: 1x2 matrices giving the indices of the hand acceleration, impact, and movement periods resp.
%       pass empty matrices to use the automatic segmentation algorithm
%   - accel, friction, sound: (structures with fields 'data' and 'Fs') the actual data
%       accel and friction .data are Nx1 arrays, sound.data is Nx2
%       instead of these three parameters, you can pass a structure with 'acc', 'fric', 'snd' fields
function [features, hac, imp, mov, timings] = steinbach_features(feats, hac, imp, mov, accel, friction, sound)
    if nargin == 5
        strukt = accel;
        accel = strukt.acc;
        friction = strukt.fric;
        sound = strukt.snd;
    end

    features = [];
    timings = containers.Map;

    %% detect impact and movement phases
    flat = [round(accel.Fs*.25) round(accel.Fs*.75)]; % assume no motion during this time
    if isempty(hac) && isempty(imp) && isempty(mov)
        flat_power = mean(abs(accel.data(flat(1):flat(2))));
        imp_start = find(abs(accel.data(flat(2):round(0.25*length(accel.data)))) > flat_power * 20, 1) + flat(2);
        if numel(imp_start) == 0
            % try again with shorter flat time assumption
            flat = [round(accel.Fs*.25) round(accel.Fs*.5)]; % assume no motion during this time
            flat_power = mean(abs(accel.data(flat(1):flat(2))));
            imp_start = find(abs(accel.data(flat(2):round(0.2*length(accel.data)))) > flat_power * 20, 1) + flat(2);
            if numel(imp_start) == 0
                error('No impact found in first 25% of data');
            end
        end

        hac = [flat(2) imp_start-1]; % TODO find start of hand acceleration
        imp = [imp_start imp_start+round(0.2*accel.Fs)];
        mov = [imp(2)+1 length(accel.data)]; % TODO find start of hand motion
    end

    %% MFCC (accel feature) -> MF
    if any(strcmp('MF', feats))
        tic;
        mf = mfcc(accel.data(mov(1):mov(2)), accel.Fs, ...
                  25, 15, ... % 25ms windows with 10ms overlap
                  0, ... % set alpha=0 to effectively remove the filtering step
                  @hamming, ... % personal correspondence
                  [0 1400], ... % personal correspondence
                  20, 20, ... % personal correspondence
                  0); % set L=0 to disable liftering

        MF = mean(mf(2:14,:)');
        features = add(features, MF);
        timings('MF') = toc;
    end
    
    %% hardness (accel feature) -> H
    % split acceleration trace into hand acceleration $h$ and impact data $i$
    % find three largest contact impulses in $i$ and calculate temporal centroid $n$
    % H = max(i)/n * 1/sum(h)
    if any(strcmp('H', feats))
        tic;
        dbclear if naninf;
        [~,n] = findpeaks(accel.data(imp(1):imp(2)), 'NPeaks',3, 'MinPeakProminence',0.01, 'MinPeakDistance',100);
        dbstop if naninf;
        acc_hand = abs(accel.data(hac(1):hac(2)));
        Hd = design(fdesign.lowpass('N,F3db', 1, 12*2*pi/accel.Fs), 'butter');
        acc_hand = filter(Hd, acc_hand);
        H = max(abs(accel.data(imp(1):imp(2)))) / mean(n) / sum(acc_hand);
        features = add(features, H);
        timings('H') = toc;
    end

    %% damping (accel feature) -> SC
    % split acceleration trace into hand acceleration $h$ and impact data $i$
    % I = DCT(i, m=4096)
    % SC = sum(abs(I(1:m/2)).^2 .* (1:m/2)) / sum(abs(I(1:m/2)).^2)
    if any(strcmp('SC', feats))
        tic;
        i = accel.data(imp(1):imp(2));
        m = 4096;
        I = dct(i, m);
        SC = sum(I(1:m/2).^2 .* linspace(0, accel.Fs/2, m/2)') / sum(I(1:m/2).^2);
        features = add(features, SC);
        timings('SC') = toc;
    end

    %% roughness (accel feature) -> TR, SR
    % high pass filter cutoff 100
    % Coiflet3 wavelet transform, extract detail levels d1 and d5
    % TR = log10(mean(d1 - mean(d1)/mean(d5)*d5))
    % take moving windows of 5000 points $x_n$, every 100 samples
    % X_n = abs(DCT(x_n))
    % D_k = mean(X_n - X_{n+100}) for each window k
    % SR = log10(sum(D_k.^2))
    if any(strcmp('TR', feats))
        tic;
        Hd = design(fdesign.highpass('N,F3db', 1, (accel.Fs/100)*2*pi/accel.Fs), 'butter');
        accfilt = filter(Hd, accel.data(mov(1):mov(2)));
        [C,L] = wavedec(accfilt, 5, 'coif3');
        d1 = abs(wrcoef('d', C, L, 'coif3', 1));
        d5 = abs(wrcoef('d', C, L, 'coif3', 5));

        avg = abs(mean(d1 - mean(d1)/mean(d5)*d5));
        if avg == 0
            TR = 0;
        else
            TR = log10(avg);
        end
        features = add(features, TR, false);
        timings('TR') = toc;
    end
    if any(strcmp('SR', feats))
        tic;
        Dksum = 0;
        shift = accel.Fs/100;
        frame = min(diff(mov)-shift, accel.Fs/2); % 1000?
        dclen = frame; % same as frame length? paper seems to say yes but personal correspondence sets frame=1000
        for n=mov(1):shift:(mov(2)-shift-frame)
            x1 = accel.data(n:n+frame);
            x2 = accel.data(n+shift:n+shift+frame);
            X1 = abs(dct(x1, dclen));
            X2 = abs(dct(x2, dclen));
            Dksum = Dksum + mean(abs(X1 - X2).^2);
        end

        SR = log10(Dksum);
        features = add(features, SR);
        timings('SR') = toc;
    end

    %% macroscopic roughness -> WV, SP, F, RG
    if any(strcmp('WV', feats))
        tic;
        % x = movement data
        % xfilt = x LPF 100 Hz
        % m = mean absolute value of each 200-sample frame of xfilt
        % s = 100-sample moving average of m (despite what the paper says)
        % WV = log10(std(m - s))
        x = accel.data(mov(1):mov(2));
        Hd = design(fdesign.lowpass('N,F3db', 1, (accel.Fs/100)*2*pi/accel.Fs), 'butter');
        xfilt = filter(Hd, x);
        framelen = accel.Fs/50;
        m = zeros(1, floor(length(xfilt)/framelen));
        for i=1:length(m)
            m(i) = mean(abs(xfilt(((i-1)*framelen+1):(i*framelen))));
        end
        s = filter(ones(100,1)/100, 1, abs(m));
        WV = 1 + log10(std(m - s)); % the +1 comes from personal correspondence, but why???
        features = add(features, WV);
        timings('WV') = toc;
    end
    if any(strcmp('SP', feats))
        tic;
        % use xfilt from WV
        % x5000 = 5000-sample moving average of xfilt
        % xth = 2*std(xfilt) + mean(xfilt) + mean(x5000)
        % xdel = max(0, xfilt - xth)
        % SP = log10(mean(xdel))
        x = accel.data(mov(1):mov(2));
        Hd = design(fdesign.lowpass('N,F3db', 1, (accel.Fs/100)*2*pi/accel.Fs), 'butter');
        xfilt = filter(Hd, x);
        x5000 = filter(ones(accel.Fs/2,1)/(accel.Fs/2), 1, xfilt);
        x5000(isnan(x5000)) = 0;
        xth = 2*std(xfilt) + mean(xfilt) + mean(x5000);
        xdel = xfilt - xth;
        xdel(xdel < 0) = 0;
        if mean(xdel) <= 0
            SP = 0;
        else
            SP = log10(mean(xdel));
        end
        features = add(features, SP, false);
        timings('SP') = toc;
    end
    if any(strcmp('F', feats))
        tic;
        % F = same as SC, but on movement data
        % the paper talks about using 1s chunks, but doesn't say how to turn that into a single feature
        % my plan is to average the F over all 1s chunks
        x = accel.data(mov(1):mov(2));
        m = 4096;
        c = accel.Fs;
        Fsum = 0;
        for i=1:c:max(1, length(x)-c)
            chunk = x(i:min(end,i+c));
            Chunk = dct(chunk, m);
            Fsum = Fsum + sum(abs(Chunk(1:m/2)).^2 .* linspace(0, accel.Fs/2, m/2)') / sum(abs(Chunk(1:m/2)).^2);
        end
        F = Fsum/length(1:c:max(1, length(x)-c));
        features = add(features, F);
        timings('F') = toc;
    end
    if any(strcmp('RG', feats))
        tic;
        x = accel.data(mov(1):mov(2));
        xhat = x/(max(x) - min(x));
        r = xcorr(xhat);
        dr = diff(r);
        dr(dr < 0) = 0;
        RG = mean(dr);
        features = add(features, RG);
        timings('RG') = toc;
    end

    %% friction features
    if size(friction.data,1) == 99991
        % scale mov/[1 acc.len] to fmov/[1 fric.len]
        flen = 99991;
        alen = size(accel.data,1) * friction.Fs/accel.Fs;
        fmov = round(flen/alen * mov);
    else
        fmov = round(mov*friction.Fs/accel.Fs);
    end
    if any(strcmp('Fr', feats)) % only in conference paper
        tic;
        Fr = mean(friction.data(fmov(1):fmov(2)));
        features = add(features, Fr);
        timings('Fr') = toc;
    end
    if any(strcmp('FM', feats)) % only in journal paper
        tic;
        vd = abs(friction.data(fmov(1):fmov(2)));
        x = accel.data(mov(1):mov(2));
        a = 1;
        b = 1;
        FM = (a*mean(vd) + b*std(diff(vd)))/mean(abs(x));
        features = add(features, FM);
        timings('FM') = toc;
    end

    %% sound features
    if any(strcmp('SIH', feats))
        tic;
        dbclear if naninf;
        [~,n] = findpeaks(accel.data(imp(1):imp(2)), 'NPeaks',3, 'MinPeakProminence',0.01, 'MinPeakDistance',100);
        dbstop if naninf;
        acc_hand = abs(accel.data(hac(1):hac(2)));
        Hd = design(fdesign.lowpass('N,F3db', 1, 12*2*pi/accel.Fs), 'butter');
        acc_hand = filter(Hd, acc_hand);
        simp = round(imp*sound.Fs/accel.Fs);
        H = max(abs(sound.data(simp(1):simp(2),1))) / mean(n) / sum(acc_hand);
        features = add(features, H);
        timings('SIH') = toc;
    end
    if any(strcmp('SMLH', feats))
        tic;
        L = size(sound.data,1);
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        smov = round(mov*sound.Fs/accel.Fs);
        SM = fft(sound.data(smov(1):smov(2),1),NFFT)/L; % FIXME why is this divided by L
        f = sound.Fs/2*linspace(0,1,NFFT/2+1);
        lo = [1 4000];
        hi = [4000 10000];
        SM_lo = SM(f >= lo(1) & f <= lo(2));
        SM_hi = SM(f >= hi(1) & f <= hi(2));
        LH = abs(mean(SM_lo))/abs(mean(SM_hi));
        features = add(features, LH);
        timings('SMLH') = toc;
    end
    if any(strcmp('SISR', feats))
        tic;
        P = pwelch(sound.data(:,1), [], [], [], sound.Fs);
        Psum = cumsum(P);
        SR = Psum(find(Psum >= 0.95*Psum(end), 1));
        features = add(features, SR);
        timings('SISR') = toc;
    end
    if any(strcmp('SISH', feats))
        tic;
        simp = round(imp*sound.Fs/accel.Fs);
        Xi = xcorr(dct(sound.data(simp(1):simp(2),1)));
        dXi = diff(Xi);
        dXi(dXi < 0) = 0;
        i = ceil(length(dXi)/2);
        f = linspace(1, sound.Fs/2, length(dXi)-i+1)';
        SH = sum(dXi(i:end) .* f)/sum(dXi(i:end));
        features = add(features, SH);
        timings('SISH') = toc;
    end
    if any(strcmp('SISS', feats)) || any(strcmp('SISC', feats))
        tic;
        % first calculate SISC
        simp = round(imp*sound.Fs/accel.Fs);
        i = sound.data(simp(1):simp(2),1);
        m = 4096;
        I = dct(i, m);
        f = linspace(0, sound.Fs/2, m/2)';
        SC = sum(I(1:m/2).^2 .* f) / sum(I(1:m/2).^2);
        % then use that to calculate SISS
        SS = sum((f - SC).^2 .* I(1:m/2).^2) / sum(I(1:m/2).^2);
        if any(strcmp('SISC', feats))
            features = add(features, SC);
        end
        if any(strcmp('SISS', feats))
            features = add(features, SS);
        end
        timings('SISS') = toc;
    end
end

function features = add(old, new, check)
    if nargin == 2 || check == true
        if isnan(new)
            fprintf('\nFeature is NaN!\n');
            keyboard;
        elseif isinf(new)
            fprintf('\nFeature is Inf!\n');
            keyboard;
        elseif new == 0
            fprintf('\nFeature is zero!\n');
            keyboard;
        end
    end

    features = [old new];
end
