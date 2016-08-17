function offset = findOffset(v,f)
%Finds the offset between the force data and the vicon data
%f is the force data, v is the vicon data
%
%If it cannot find a reasonable offset, then the offset will be set to 0
%      (Find a better number for this)

    if isempty(v)
        offset = 0;
        return;
    end
    
    % find the 4 taps in force (hand-tuned parameters seem to work)
    [~, flocs] = findpeaks(f(:,3), 'SortStr','descend', 'NPeaks',4, 'MinPeakProminence',2, 'MinPeakDistance',500, 'MinPeakHeight',27);
    flocs = sort(flocs); % resort in order
    ftimes = f(flocs([1 3]),1); % use the first and second-to-last
    
    % find the 4 rises-before-the-taps in vicon
    [~, vlocs] = findpeaks(v(:,4), 'SortStr','descend', 'NPeaks',4, 'MinPeakDistance',10);
    vlocs = sort(vlocs); % resort in order
    [~,i1] = min(v(vlocs(1):vlocs(2),4));
    [~,i2] = min(v(vlocs(3):vlocs(4),4));
    vtimes = [v(i1 + vlocs(1) - 1, 1); v(i2 + vlocs(3) - 1, 1)]; % find min between each pair to get the first and second-to-last
    
    % the offset is the average time difference between the peaks,
    % but we might have some misidentified peaks.
    offset = mean(vtimes - ftimes);
    
    if abs(offset) > 20
        offset = 0;
        disp('Could not find appropriate offset. Data may not be aligned.')
    end
end