function [v, int, dig_acc, dig_gyro, mic, ana_acc, mag, dt, opto, bio] = load_stick(prefix)
%v, vicon (x, y, z, rotation)
%int force (internal midi40)
%dig_acc and dig_gyro: IMU data
%mic: microphone
%ana_acc: accelerometers
%mag: magnatometer (IMU)
%raw: raw force data (not going to be using this, most likely)
    v = csvload([prefix 'vicon.tsv'], ...
                {'Timestamp', ...
                 'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
                 'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
                {'Delimiter', '\t'});

    if nargout > 1
        try
            int = csvload([prefix 'teensy.ft.csv'], ...
                          [{'TeensyDt', 'Timestamp'}, ...
                           arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
            dt = int(:,1);
            int = int(:,2:end);
        catch err
            if strcmp(err.message, 'Could not find column TeensyDt')
                % old format
                int = csvload([prefix 'teensy.ft.csv'], ...
                              [{'Timestamp'}, ...
                              arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
            end
        end

        acc = csvload([prefix 'teensy.acc.csv'], ...
                      {'Timestamp', 'FIFOPosition', 'AccX', 'AccY', 'AccZ'});
        gyro = csvload([prefix 'teensy.gyro.csv'], ...
                      {'Timestamp', 'FIFOPosition', 'GyroX', 'GyroY', 'GyroZ'});
        mag = csvload([prefix 'teensy.mag.csv'], ...
                      {'Timestamp', 'MagX', 'MagY', 'MagZ'});

        % unpack raw sensor data
        [int, mic, ana_acc] = process_mini40(int);%, zeros(1,6), eye(6,6));
        dig_acc = unfifo(acc);
        dig_gyro = unfifo(gyro);
        
        % load optoforce and biotac if present
        opto = csvload([prefix 'optoforce.csv'], ...
                       {'Timestamp', 'X', 'Y', 'Z'});
        bio = csvload([prefix 'biotac.csv'], ...
                       {'Timestamp', 'PDC', 'PAC_0', 'PAC_1', 'PAC_2', 'PAC_3', 'PAC_4', 'PAC_5', 'PAC_6', 'PAC_7', 'PAC_8', 'PAC_9', 'PAC_10', 'PAC_11', 'PAC_12', 'PAC_13', 'PAC_14', 'PAC_15', 'PAC_16', 'PAC_17', 'PAC_18', 'PAC_19', 'PAC_20', 'PAC_21', 'TDC', 'TAC', 'Electrode_0', 'Electrode_1', 'Electrode_2', 'Electrode_3', 'Electrode_4', 'Electrode_5', 'Electrode_6', 'Electrode_7', 'Electrode_8', 'Electrode_9', 'Electrode_10', 'Electrode_11', 'Electrode_12', 'Electrode_13', 'Electrode_14', 'Electrode_15', 'Electrode_16', 'Electrode_17', 'Electrode_18'});
    end
end
