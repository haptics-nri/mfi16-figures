function [start, stop, spike_clumps] = narrow_to_taps(force)    

    % start 0.5s after last beginning spike, end 0.5s before first ending spike
    spikes = find(diff(force(:,4)) > 0.5);
    spike_boundaries = [1 find(diff(spikes) > 100)' length(spikes)];
    spike_clumps = zeros([1 length(spike_boundaries)-1]);
    for i=1:length(spike_clumps)
        spike_clumps(i) = mean(spikes(spike_boundaries(i):spike_boundaries(i+1)));
    end
    early_taps = [1 spike_clumps(spike_clumps < size(force,1)*1/3)]; % FIXME this is fragile
    late_taps = [spike_clumps(spike_clumps > size(force,1)*2/3) size(force,1)];
    start = round(early_taps(end) + 3000);
    stop = round(late_taps(1) - 1500);
    
end
