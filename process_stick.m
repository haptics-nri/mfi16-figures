function [v, int, vbody, vend, vint, vbodyint, vendint, accint, accworld, intbody, intworld, intworldsub] = process_stick(v, int, acc, mass, com, Hvb, Hmb, Hib, off)

    % narrow to common time window
    fprintf('\tnarrowing\n');
    start = max([int(1,1)-off, v(1,1)]);
    stop = min([int(end,1)-off, v(end,1)]);
    int = int(int(:,1)-off>=start+.03 & int(:,1)-off<=stop-.03, :); %%% HACK
    v = v(v(:,1)>=start & v(:,1)<=stop, :);
    
    % remove spikes (HACK)
    fprintf('\tsmoothing timestamps\n');
    spikes = find(diff(int(:,1)) > 2*mean(diff(int(:,1))));
    for i=1:length(spikes)
        N = 10;
        int(spikes(i)-N:spikes(i)+N,1) = int(spikes(i)-N) + ((int(spikes(i)+N) - int(spikes(i)-N)) * linspace(0,1,2*N+1));
    end

    % transform Vicon into body frame and end-effector frame
    fprintf('\ttransforming pose\n');
    vbody = v;
    vend = v;
    for i=1:size(vbody,1)
        tf = xfconv(v(i,2:7)) * Hvb;
        vbody(i,2:4) = tf(1:3,4);
        vbody(i,5:7) = xfconv(tf(1:3,1:3));
        tf = tf * Hib;
        vend(i,2:4) = tf(1:3,4);
        vend(i,5:7) = xfconv(tf(1:3,1:3));
    end

    % upsample and offset Vicon/accel to match Mini40
    fprintf('\tupsampling\n');
    vint     = [int(:,1) interp1(v(:,1),     v(:,2:4),     int(:,1)-off) slerp(v(:,1),     v(:,5:7),     int(:,1)-off)];
    vbodyint = [int(:,1) interp1(vbody(:,1), vbody(:,2:4), int(:,1)-off) slerp(vbody(:,1), vbody(:,5:7), int(:,1)-off)];
    vendint  = [int(:,1) interp1(vend(:,1),  vend(:,2:4),  int(:,1)-off) slerp(vend(:,1),  vend(:,5:7),  int(:,1)-off)];
    accint   = [int(:,1) interp1(acc(:,1),   acc(:,2:4),   int(:,1))];

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
        accworld(i,2:4) = tf * accint(i,2:4)';
    end

    % subtract weight of end-effector
    fprintf('\tgravity compensation\n');
    intworldsub = intbody;
    for i=1:size(intbody,1)
        fg = mass * [0; 0; -9.81];
        intworldsub(i,2:4) = intworld(i,2:4) - fg';
        intworldsub(i,5:7) = intworld(i,5:7) - cross(com, fg)';
    end
    
    % HACK smooth complex parts of intworldsub
    

end
