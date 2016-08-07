% part 4 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    j = str2double(epdirs(i).name);
    prefix = [datadir filesep epdirs(i).dir filesep 'stick' filesep epdirs(i).name];
    fprintf('[%d/%d] %s\n', i, length(epdirs), prefix);
    eps(j).flow = parse_flow([prefix filesep 'stick.flow']);
    eps(j).endeff = eps(j).flow.answers('tooling ball diameter').text;
    eps(j).material = eps(j).flow.answers('surface name').text;
    [v, f, da, dg, ~, a] = load_stick([prefix filesep]);
    eps(j).data = struct('vicon', v, 'force', f, 'acc', a, 'imu', struct('acc', da, 'gyro', dg));
