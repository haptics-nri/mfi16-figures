function data = unfifo(data)

    i = 1;
    while i < size(data,1)
        % j is the start of the next packet
        j = i + 1;
        while j < size(data,1) && data(j,2) ~= 0
            j = j + 1;
        end
        
        % interpolate timestamps through the packet
        data(i:j,1) = linspace(data(i,1), data(j,1), j-i+1);
        
        % jump to next packet (record length)
        i = j;
    end
    
    % remove rows with repeated time values (otherwise interp1 fails)
    data = data(diff(data(:,1)) > 0, :);
    
    % remove FIFOPosition column
    data = data(:, [1 3:end]);

end
