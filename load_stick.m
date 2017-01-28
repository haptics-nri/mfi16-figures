function [v, int, dig_acc, dig_gyro, mic, ana_acc, mag, dt, opto, bio, motrak] = load_stick(prefix)
%v, vicon (x, y, z, rotation)
%int force (internal midi40)
%dig_acc and dig_gyro: IMU data
%mic: microphone
%ana_acc: accelerometers
%mag: magnatometer (IMU)
%raw: raw force data (not going to be using this, most likely)

    if prefix(end) ~= '/'
        prefix = [prefix '/'];
    end

    % load vicon (if present)
    v = csvload([prefix 'vicon.tsv'], ...
                {'Timestamp', ...
                 'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
                 'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
                {'Delimiter', '\t'});

    if nargout > 1
        ftname = 'teensy.ft.csv';
        if exist([prefix 'teensy.ft.fixed.csv'], 'file')
            ftname = 'teensy.ft.fixed.csv';
        end
        
        try
            int = csvload([prefix ftname], ...
                          [{'TeensyDt', 'PacketNumber', 'Timestamp'}, ...
                           arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
            dt = int(:,1);
            accref = int(:,2);
            int = int(:,3:end);
        catch err
            if strcmp(err.message, 'Could not find column TeensyDt')
                % old format
                int = csvload([prefix 'teensy.ft.csv'], ...
                              [{'Timestamp'}, ...
                              arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
            end
            dt = [];
            accref = 100*ones(size(int,1),1);
        end

        acc = csvload([prefix 'teensy.acc.csv'], ...
                      {'Timestamp', 'FIFOPosition', 'AccX', 'AccY', 'AccZ'});
        gyro = csvload([prefix 'teensy.gyro.csv'], ...
                      {'Timestamp', 'FIFOPosition', 'GyroX', 'GyroY', 'GyroZ'});
        mag = csvload([prefix 'teensy.mag.csv'], ...
                      {'Timestamp', 'MagX', 'MagY', 'MagZ'});

        % unpack raw sensor data
        [int, mic, ana_acc] = process_mini40(accref, int);%, zeros(1,6), eye(6,6));
        dig_acc = unfifo(acc);
        dig_gyro = unfifo(gyro);
        
        % load optoforce and biotac if present
        opto = csvload([prefix 'optoforce.csv'], ...
                       {'Timestamp', 'X', 'Y', 'Z'});
        bio = csvload([prefix 'biotac.csv'], ...
                       {'Timestamp', 'PDC', 'PAC_0', 'PAC_1', 'PAC_2', 'PAC_3', 'PAC_4', ...
                        'PAC_5', 'PAC_6', 'PAC_7', 'PAC_8', 'PAC_9', 'PAC_10', 'PAC_11', ...
                        'PAC_12', 'PAC_13', 'PAC_14', 'PAC_15', 'PAC_16', 'PAC_17', 'PAC_18', ...
                        'PAC_19', 'PAC_20', 'PAC_21', 'TDC', 'TAC', 'Electrode_0', ...
                        'Electrode_1', 'Electrode_2', 'Electrode_3', 'Electrode_4', ...
                        'Electrode_5', 'Electrode_6', 'Electrode_7', 'Electrode_8', ...
                        'Electrode_9', 'Electrode_10', 'Electrode_11', 'Electrode_12', ...
                        'Electrode_13', 'Electrode_14', 'Electrode_15', 'Electrode_16', ...
                        'Electrode_17', 'Electrode_18'});
        
        % load april tag data from bluefox (if present)
        ts_rgb = csvload([prefix 'bluefox/bluefox_times.csv'], {'UnixTimestamp'});
        motrak = [];
        if ~isempty(ts_rgb)
            % if there's april tag data, do motion tracking on it

            [nums, ids, centers, p1s, p2s, p3s, p4s] = load_april([prefix, 'bluefox/', 'april.csv']);
            [~,order] = sort(nums);
            centers = centers(order);
            ids = ids(order);
            p1s = p1s(order);
            p2s = p2s(order);
            p3s = p3s(order);
            p4s = p4s(order);

            pos = zeros(length(nums),3);
            motrak = [ts_rgb zeros(length(nums),6)];

            for i = 1:length(nums)
                ts = ts_rgb(i);
                data(i).t = ts;
                data(i).id = ids{i};
                data(i).p0 = centers{i};
                data(i).p1 = p1s{i};
                data(i).p2 = p2s{i};
                data(i).p3 = p3s{i};
                data(i).p4 = p4s{i};

                %%% Uncomment it for fusion of PNP + IMU
                %[~,idx_omg] = min(abs(gyro(:,1)-ts));
                %data(i).omg = (omegas(idx_omg,:))';
                %[~,idx_acc] = min(abs(acc(:,1)-ts));
                %data(i).acc = (accelerations(idx_acc,:))';

                [X,Q] = estimate_pose_pnp(data(i));

                %%% Uncomment it for fusion of PNP + IMU
                %[X,~] = NewFusion(data(i));
                %Q = X(7:10);

                %%% Calculating position of end effector from the position of IMU
                cal_results;
                R_body2imu =  [1, 0,  0;
                              0,-1,  0;
                              0, 0, -1];

                T_imu2 = [-254.3402 + 191.08152942; 0; 30.96486407];


                H_body2imu = [R_body2imu,T_imu2;0,0,0,1];
                H_matrix = [quat2rotm(Q'),X(1:3);[0,0,0,1]];
                newMatrix = H_matrix*H_body2imu*H_bal2imu;

                pos(i,:) = newMatrix(1:3,4);
                
                motrak(i,2:4) = pos(i,:);
                motrak(i,5:7) = xfconv(newMatrix(1:3,1:3));
            end


            % if there's no vicon but there is bluefox, let's make believe
            if isempty(v) && ~isempty(motrak)
                v = motrak;
            end
            
        end
        
    end
end
