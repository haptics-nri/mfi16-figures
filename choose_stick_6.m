% part 6 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    fprintf('%d\n', i);
    
    v = eps(i).data.vicon;
    f = eps(i).data.force;
    
    % find the 4 taps in force (hand-tuned parameters seem to work)
    [fpks, flocs] = findpeaks(f(:,3), 'SortStr','descend', 'NPeaks',4, 'MinPeakProminence',2, 'MinPeakDistance',500);
    [flocs, idx] = sort(flocs); % resort in order
    fpks = fpks(idx);
    ftimes = f(flocs([1 3]),1); % use the first and second-to-last
    
    % find the 4 rises-before-the-taps in vicon
    [vpks, vlocs] = findpeaks(v(:,4), 'SortStr','descend', 'NPeaks',4, 'MinPeakDistance',10);
    [vlocs, idx] = sort(vlocs); % resort in order
    vpks = vpks(idx);
    [~,i1] = min(v(vlocs(1):vlocs(2),4));
    [~,i2] = min(v(vlocs(3):vlocs(4),4));
    vtimes = [v(i1 + vlocs(1) - 1, 1); v(i2 + vlocs(3) - 1, 1)]; % find min between each pair to get the first and second-to-last
    
    % the offset is the average time difference between the peaks,
    % but we might have some misidentified peaks.
    eps(i).offset = mean(vtimes - ftimes);
    eps(i).peaks = flocs;
