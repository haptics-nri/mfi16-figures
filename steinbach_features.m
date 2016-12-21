function [features, imp, mov] = steinbach_features(feats, accel, friction, sound)
    if nargin == 2
        strukt = accel;
        accel = strukt.acc;
        friction = strukt.fric;
        sound = strukt.snd;
    end

    features = [];

    %% detect impact and movement phases
    flat = [round(accel.Fs*.25) round(accel.Fs*.75)]; % assume no motion during this time
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
    imp = [imp_start imp_start+round(0.2*accel.Fs)];
    mov = [imp(2)+1 length(accel.data)]; % TODO find start of hand motion

    %% MFCC (accel feature) -> MF
    if any(strcmp('MF', feats))
        mf = mfcc(accel.data, accel.Fs, ...
                  25, 15, ... % 25ms windows with 10ms overlap
                  0, ... % set alpha=0 to effectively remove the filtering step
                  @hamming, ... % personal correspondence
                  [0 1400], ... % personal correspondence
                  20, 20, ... % personal correspondence
                  0); % set L=0 to disable liftering

        MF = mean(mf(2:14,:)');
        features = [features MF];
    end
    
    %% hardness (accel feature) -> H
    % split acceleration trace into hand acceleration $h$ and impact data $i$
    % find three largest contact impulses in $i$ and calculate temporal centroid $n$
    % H = max(i)/n * 1/sum(h)
    % but how is the spike detection done???

    %% damping (accel feature) -> SC
    % split acceleration trace into hand acceleration $h$ and impact data $i$
    % I = DCT(i, m=4096)
    % SC = sum(abs(I(1:m/2)).^2 .* (1:m/2)) / sum(abs(I(1:m/2)).^2)
    if any(strcmp('SC', feats))
        i = accel.data(imp(1):imp(2));
        m = 4096;
        I = dct(i, m);
        SC = sum(abs(I(1:m/2)).^2 .* linspace(0, accel.Fs/2, m/2)') / sum(abs(I(1:m/2)).^2);
        features = [features SC];
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
        Hd = design(fdesign.highpass('N,F3db', 1, 300*2*pi/accel.Fs), 'butter');
        accfilt = filter(Hd, accel.data(mov(1):mov(2)));
        [C,L] = wavedec(accfilt, 5, 'coif3');
        d1 = abs(wrcoef('d', C, L, 'coif3', 1));
        d5 = abs(wrcoef('d', C, L, 'coif3', 5));

        TR = log10(abs(mean(d1 - mean(d1)/mean(d5)*d5)));
        features = [features TR];
    end
    if any(strcmp('SR', feats))
        Dksum = 0;
        shift = 100;
        frame = 5000; % 1000?
        dclen = 5000; % same as frame length? paper seems to say yes but personal correspondence sets frame=1000
        for n=mov(1):shift:(mov(2)-shift-frame)
            x1 = accel.data(n:n+frame);
            x2 = accel.data(n+shift:n+shift+frame);
            X1 = abs(dct(x1, dclen));
            X2 = abs(dct(x2, dclen));
            Dksum = Dksum + mean(abs(X1 - X2).^2);
        end

        SR = log10(Dksum);
        features = [features SR];
    end

    %% macroscopic roughness -> WV, SP, F, RG
    if any(strcmp('WV', feats))
        % x = movement data
        % xfilt = x LPF 100 Hz
        % m = mean absolute value of each 200-sample frame of xfilt
        % s = 100-sample moving average of m (despite what the paper says)
        % WV = log10(std(m - s))
        x = accel.data(mov(1):mov(2));
        Hd = design(fdesign.lowpass('N,F3db', 1, 100*2*pi/accel.Fs), 'butter');
        xfilt = filter(Hd, x);
        m = zeros(1, floor(length(xfilt)/200));
        for i=1:length(m)
            m = mean(abs(xfilt(((i-1)*200+1):(i*200))));
        end
        s = filter(ones(100,1)/100, 1, abs(m));
        WV = 1 + log10(std(m - s)); % the +1 comes from personal correspondence, but why‽‽‽
        features = [features WV];
    end
    if any(strcmp('SP', feats))
        % use xfilt from WV
        % x5000 = 5000-sample moving average of xfilt
        % xth = 2*std(xfilt) + mean(xfilt) + mean(x5000)
        % xdel = max(0, xfilt - xth)
        % SP = log10(mean(xdel))
        x = accel.data(mov(1):mov(2));
        Hd = design(fdesign.lowpass('N,F3db', 1, 100*2*pi/accel.Fs), 'butter');
        xfilt = filter(Hd, x);
        x5000 = tsmovavg(xfilt, 's', 5000, 1);
        x5000(isnan(x5000)) = 0;
        xth = 2*std(xfilt) + mean(xfilt) + mean(x5000);
        xdel = xfilt - xth;
        xdel(xdel < 0) = 0;
        SP = log10(mean(xdel));
        features = [features SP];
    end
    if any(strcmp('F', feats))
        % F = same as SC, but on movement data
        % the paper talks about using 1s chunks, but doesn't say how to turn that into a single feature
        % my plan is to average the F over all 1s chunks
        x = accel.data(mov(1):mov(2));
        m = 4096;
        c = accel.Fs;
        Fsum = 0;
        for i=1:c:(length(x)-c)
            chunk = x(i:i+c);
            Chunk = dct(chunk, m);
            Fsum = Fsum + sum(abs(Chunk(1:m/2)).^2 .* linspace(0, accel.Fs/2, m/2)') / sum(abs(Chunk(1:m/2)).^2);
        end
        F = Fsum/length(1:c:length(x)-c);
        features = [features F];
    end
    if any(strcmp('RG', feats))
        x = accel.data(mov(1):mov(2));
        xhat = x/(max(x) - min(x));
        r = xcorr(xhat);
        dr = diff(r);
        dr(dr < 0) = 0;
        RG = mean(dr);
        features = [features RG];
    end

    %% sound features
end

