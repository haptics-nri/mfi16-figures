function ft = process_mini40(raw, bias, tf)

    if nargin < 3
        % from mfgr on proton pack internal Mini40
        tf =           [ 0.00679   0.01658  -0.04923   6.20566   0.15882  -6.19201 
                         0.11638  -7.31729  -0.04322   3.54949  -0.08024   3.57115 
                        10.35231   0.32653  10.61091   0.29668  10.33382   0.25761 
                         0.00022  -0.04140   0.14917   0.02435  -0.15234   0.01567 
                        -0.16837  -0.00464   0.08561  -0.03311   0.08763   0.03721 
                         0.00128  -0.08962   0.00085  -0.08785   0.00204  -0.08790 ];
        if nargin < 2
            % measured on proton pack internal Mini40
            x = load('mini40_calib.mat');
            bias = (tf\[x.fbias x.tbias]')';
        end
    end
    
    % merge bytes (preserve timestamp column), handle negatives, and scale
    ft = [raw(:,1) bitshift(raw(:,[2 4 6 8 10 12]), 8) + raw(:,[3 5 7 9 11 13])];
    for i=1:size(ft,1)
        for j=2:size(ft,2)
            if ft(i,j) >= 2048
                ft(i,j) = ft(i,j) - 4096;
            end
        end
    end
    ft(:,2:end) = ft(:,2:end) * 0.002;
    ft = permute(ft, [1 7 6 5 4 3 2]); % why is this necessary?

    % bias and transform
    ft(:,2:end) = (tf * bsxfun(@minus, ft(:,2:end), bias)')';

    ft = squeeze(ft);
end
