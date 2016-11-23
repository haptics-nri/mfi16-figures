function [v, int, vbody, vend, vint, vbodyint, vendint, accint, accworld, intbody, intworld, intworldsub, start, stop] = process_stick(v, int, acc, mass, com, Hvb, Hmb, Hib, off)
% Inputs:
% v: vicon data (nx7; first column time)
% int: force (time, force x, y, z torque x, y, z)
% acc: vibrations (nx7: time, x, y, z from first, x, y, z from second)
% mass/com/transformation matrices
% off: time offset
% 
% Outputs:
% v: Position of marker frame (time has been narrowed)
% int: force
% vbody: position of the body frame in world frame
% vend: position of end in the world frame
% ( )-int: interpolated so that sampling rate matches sampling rate of the
%   force
% accworld: acc in the world frame *plot this*
% intbody: force in body frame
% intworld: force in the world frame
% intworldsub: force in the world frame compensated for gravity


    if 1/mean(diff(v(:,1))) < 50
        is_bluefox = true;
    else
        is_bluefox = false;
    end

    if nargin == 8 || isempty(off)
        if is_bluefox
            off = 0;
        else
            %Testing out this offset thing
            off = -findOffset(v, int, 100)
        end
    end

    if ~isempty(v)
        % narrow to common time window
        fprintf('\tnarrowing\n');
        start = max([int(1,1)-off, v(1,1)]);
        stop = min([int(end,1)-off, v(end,1)]);
        int = int(int(:,1)-off>=start & int(:,1)-off<=stop, :);
        v = v(v(:,1)>=start & v(:,1)<=stop, :);
    end
    
    % remove spikes (HACK)
    fprintf('\tsmoothing timestamps\n');
    spikes = find(diff(int(:,1)) > 2*mean(diff(int(:,1))));
    for i=1:length(spikes)
        % average up to 10 values from the left and right (but don't overrun the array)
        Nl = min([10 spikes(i)-1]);
        Nr = min([10 size(int,1)-spikes(i)]);
        int(spikes(i)-Nl:spikes(i)+Nr,1) = int(spikes(i)-Nl) + ((int(spikes(i)+Nr) - int(spikes(i)-Nl)) * linspace(0,1,Nl+Nr+1));
    end
    
    % smooth Vicon angles
    if ~is_bluefox
        fprintf('\tsmoothing vicon angles\n');
        v(:,5:7) = filtfilt(.4, [1 .4-1], v(:,5:7));
    end

    % transform Vicon into body frame and end-effector frame
    fprintf('\ttransforming pose\n');
    vbody = v;
    vend = v;
    if ~is_bluefox
        for i=1:size(vbody,1)
            tf = xfconv(v(i,2:7)) * Hvb;
            vbody(i,2:4) = tf(1:3,4);
            vbody(i,5:7) = xfconv(tf(1:3,1:3));
            tf = tf * Hib;
            vend(i,2:4) = tf(1:3,4);
            vend(i,5:7) = xfconv(tf(1:3,1:3));
        end
    end

    % upsample and offset Vicon/accel to match Mini40
    fprintf('\tupsampling\n');
    if isempty(v)
        vint = [];
        vbodyint = [];
        vendint = [];
    else
        vint     = [int(:,1) interp1(v(:,1),     v(:,2:4),     int(:,1)-off) slerp(v(:,1),     v(:,5:7),     int(:,1)-off)];
        vbodyint = [int(:,1) interp1(vbody(:,1), vbody(:,2:4), int(:,1)-off) slerp(vbody(:,1), vbody(:,5:7), int(:,1)-off)];
        vendint  = [int(:,1) interp1(vend(:,1),  vend(:,2:4),  int(:,1)-off) slerp(vend(:,1),  vend(:,5:7),  int(:,1)-off)];
    end
    if isempty(acc)
        accint = [];
    else
        % with the ADXL335 it must be [2 4 5] (those are the only pins that
        % are connected)
        % with the ADXL326 I am confused as to whether it should have been
        % [2 4 5] or [2 3 6] (maybe it doesn't matter?)
        accint = [int(:,1) interp1(acc(:,1), acc(:,[2 4 5]), int(:,1))];
    end
    % transform Mini40 and IMU into body frame and world frame
    fprintf('\ttransforming force\n');
    intbody = int;
    intworld = int;
    accworld = accint;
    for i=1:size(intbody,1)
        intbody(i,2:4) = Hmb(1:3,1:3)*int(i,2:4)';
        intbody(i,5:7) = Hmb(1:3,1:3)*int(i,5:7)';
        tf = xfconv(vbodyint(i,5:7));
        intworld(i,2:4) = tf * intbody(i,2:4)';
        intworld(i,5:7) = tf * intbody(i,5:7)';
        if ~isempty(accint)
            accworld(i,2:4) = tf * accint(i,2:4)';
        end
    end

    % subtract weight of end-effector
    fprintf('\tgravity compensation\n');
    intworldsub = intbody;
    for i=1:size(intbody,1)
        fg = mass * [0; 0; -9.81];
        intworldsub(i,2:4) = intworld(i,2:4) - fg';
        intworldsub(i,5:7) = intworld(i,5:7) - cross(com, fg)'; % FIXME this needs to be done in the body frame
    end
    %intworldsub(:,2:4) = filtfilt(.04, [1 .04-1], intworldsub(:,2:4));
    
end
