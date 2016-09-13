% part 9 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    if episodes(ep).name(1) == '.'
        icra17_figures_10
    elseif str2num(episodes(ep).name) < 7 % change this to select end-effector
        icra17_figures_11
    end
    
    flow = parse_flow([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep 'stickcam.flow']);
    material = flow.answers('surface name').text;
    fprintf('Loading %s (%s)...\n', episodes(ep).name, material);
    [v, int, da,dg,mi, acc, ma, dt, o,b, motrak] = ...
        load_stick([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep]);
    data(material) = struct('v',v, 'int',int, 'acc',acc, ...
                            'imu', struct('acc',da, 'gyro',dg, 'mag',ma), ...
                            'sensor', struct('opto',o, 'bio',b), ...
                            'motrak',motrak, 'dt',dt);
