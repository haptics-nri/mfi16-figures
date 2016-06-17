function [ft, mic, acc, raw] = process_mini40(raw, bias, tf)

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
    
    % COLUMNS
    % sensor        position   bytes
    % FT            1-6        2-13
    % microphone    7          14-15
    % acc 1 Y       10         20-21
    % acc 1 Z       11         22-23
    % acc 1 X       12         24-25
    % acc 2 Y       13         26-27
    % acc 2 Z       14         28-29
    % acc 2 X       15         30-31
    
    % merge bytes (preserve timestamp column), handle negatives, and scale
    ft  = [raw(:,1) bitshift(raw(:,[2 4 6 8 10 12]),     8) + raw(:,[3 5 7 9 11 13])    ];
    mic = [raw(:,1) bitshift(raw(:,14),                  8) + raw(:,15)                 ];
    acc = [raw(:,1) bitshift(raw(:,[20 22 24 26 28 30]), 8) + raw(:,[21 23 25 27 29 31])];
    ft = convsign(ft);
    ft(:,2:end) = ft(:,2:end) * 0.002;
    acc(:,2:end) = (acc(:,2:end) - 2048)/4096 * (16*9.81);

    % bias and transform
    ft(:,2:end) = (tf * bsxfun(@minus, ft(:,2:end), bias)')';
end

function cols = convsign(cols)
    for i=1:size(cols,1)
        for j=2:size(cols,2)
            if cols(i,j) >= 2048
                cols(i,j) = cols(i,j) - 4096;
            end
        end
    end
end
