function [v, int, dig_acc, dig_gyro, mic, ana_acc, mag, dt, raw] = load_stick(prefix)
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
    [int, mic, ana_acc, raw] = process_mini40(int);%, zeros(1,6), eye(6,6));
    dig_acc = unfifo(acc);
    dig_gyro = unfifo(gyro);

end
