% For use in calibrating the Proton Pack.
% Determines the mass and center-of-mass of an end-effector from
% data gathered from free-space gyration of the Proton.
% Calibrate the Mini40 bias first (and set the values in process_mini40.m).
% Add the RANSAC Toolbox to the path before calling this function.
%
% Parameters:
%    epdirs: cell array of directories containing episode files (no trailing slash)
%    display [optional, default=false]: plot something
%    testidx [optional]: pick some episodes out of epdirs to use for confirmation
%    offsets [optional, required for testidx]: clock offsets for the test episodes
% Returns:
%    mass: mass of the end-effector (kg)
%    bias: measured force with no load (N)
function [mass, bias, err] = weigh(epdirs, display, testidx, offsets)

    if nargin < 4
        if nargin == 3
            error('testidx requires offsets');
        end
        offsets = [];
        testidx = [];
        if nargin < 2
            display = false;
        end
    end

    % load data

    testint = [];
    allint = [];
    for i=1:length(epdirs)
        int = csvload([epdirs{i} filesep 'teensy.ft.csv'], ...
                      [{'Timestamp'}, ...
                       arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);

        int = process_mini40(int, zeros(1,6));
        
        a = round(size(int,1)*1/10); % middle 80% (FIXME magic number)
        b = round(size(int,1)*9/10);
        int = int(a:b,:);
        
        allint = [allint; int];
        if find(testidx == i, 1)
            testint = [testint; int];
        end
    end
    
    % perform calibration
    
    [bias, weight, ~, r_inliers] = sphereFit_ransac(allint(:,2:4));
    [~, ~, ~, r_err] = sphereFit(allint(r_inliers,2:4)); % run again to get RMS error
    mass = weight/9.81;
    err = struct('err', r_err, 'inliers', r_inliers);
    
    % display
    
    if display
        clf;
        plot(allint(:,2:4));
        sphereplot(bias, weight, {allint(r_inliers,2:4), allint(~r_inliers,2:4)});
        %t = 1:size(allint,1);
        %plot(t(r_inliers),   sqrt(sum(bsxfun(@minus, allint(r_inliers,2:4),  bias).^2,2)) - weight, '.', ...
        %     t(~r_inliers),  sqrt(sum(bsxfun(@minus, allint(~r_inliers,2:4), bias).^2,2)) - weight, '.', ...
        %     t,            acos(allint(:,3)-mean(allint(:,3))));
    end
     
    % test
    
    if ~isempty(testidx)
        testv = [];
        for i=1:length(testidx)
            v = csvload([epdirs{testidx(i)} filesep 'vicon.tsv'], ...
                        {'Timestamp', ...
                         'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
                         'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
                        {'Delimiter', '\t'});
            a = round(size(v,1)*1/10); % middle 80% (FIXME magic number)
            b = round(size(v,1)*9/10);
            v = v(a:b,:);
            v(:,1) = v(:,1) + offsets(i);
            testv = [testv; v];
        end
        v = testv;
        int = testint;
        % constants
        Hva = [ 1       0       0         5.65
                0      -1       0         6.57
                0       0      -1      -559.01
                0       0       0         1   ];
        Hib = [ 0       0       1       108.99
               -1       0       0         0.53
                0      -1       0        -2.98
                0       0       0         1   ];
        % transform Vicon into body frame and end-effector frame
        vbody = v;
        vend = v;
        for i=1:size(vbody,1)
            tf = xfconv(v(i,2:7)) * Hva;
            vbody(i,2:4) = tf(1:3,4);
            vbody(i,5:7) = xfconv(tf(1:3,1:3));
            tf = tf * Hib;
            vend(i,2:4) = tf(1:3,4);
            vend(i,5:7) = xfconv(tf(1:3,1:3));
        end

        % upsample and offset Vicon/accel to match Mini40
        vbodyint = [int(:,1) interp1(vbody(:,1), vbody(:,2:4), int(:,1)) slerp(vbody(:,1), vbody(:,5:7), int(:,1))];

        % transform Mini40 and IMU into body frame and world frame
        intbody = int;
        intworld = int;
        for i=1:size(intbody,1)
            intbody(i,2:4) = Him(1:3,1:3)*int(i,2:4)';
            intbody(i,5:7) = Him(1:3,1:3)*int(i,5:7)';
            tf = xfconv(vbodyint(i,5:7));
            intworld(i,2:4) = tf * intbody(i,2:4)';
            intworld(i,5:7) = tf * intbody(i,5:7)';
        end

        % subtract weight of end-effector
        intworldsub = intbody;
        for i=1:size(intbody,1)
            fg = mass * [0; 0; -9.81];
            intworldsub(i,2:4) = intworld(i,2:4) - fg';
            intworldsub(i,5:7) = intworld(i,5:7) - cross(com, fg)';
        end
        
        
        plot(intworldsub(:,2:4));
    end

end
