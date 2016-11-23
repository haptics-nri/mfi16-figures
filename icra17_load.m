function [data, materials] = icra17_load(DATADIR, date, flowtype, lambda)
    data = containers.Map;
    episodes = dir([DATADIR filesep date filesep flowtype]);
    for ep = 1:length(episodes)
        if episodes(ep).name(1) == '.'
            continue;
        elseif lambda(str2num(episodes(ep).name))
            continue;
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
    end

    materials = data.keys;
end
