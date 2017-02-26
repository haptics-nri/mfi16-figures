function [data, materials, calib] = icra17_load(DATADIR, dates, flowtype, lambda, weighing)
    if ischar(dates)
        dates = {dates};
    end

    if nargin == 5
        weighdate = weighing{1};
        weighep = weighing{2};
        fprintf('Weighing end-effector using %s/weigh/%s...\n', weighdate, weighep);
        [mass, fbias, ~, com, tbias] = weigh({fullfile(DATADIR, weighdate, 'weigh', weighep)});
        calib = struct('mass', mass, ...
                       'com',  com, ...
                       'fbias', fbias, ...
                       'tbias', tbias');
        varcalib = {calib};
    else
        varcalib = {};
    end

    data = containers.Map;
    for d=1:length(dates)
        date = dates{d};
        episodes = dir([DATADIR filesep date filesep flowtype]);
        for ep = 1:length(episodes)
            if episodes(ep).name(1) == '.'
                continue;
            elseif lambda(str2num(episodes(ep).name))
                continue;
            end

            flow = parse_flow([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep flowtype '.flow']);
            material = flow.answers('surface name').text;
            fprintf('Loading %s (%s)...\n', episodes(ep).name, material);
            [v, int, da,dg,mi, acc, ma, dt, o,b, motrak] = ...
                load_stick([DATADIR filesep date filesep flowtype filesep episodes(ep).name filesep], varcalib{:});
            data(material) = struct('v',v, 'int',int, 'acc',acc, ...
                                    'imu', struct('acc',da, 'gyro',dg, 'mag',ma, 'mic',mi), ...
                                    'sensor', struct('opto',o, 'bio',b), ...
                                    'motrak',motrak, 'dt',dt);
        end
    end

    materials = data.keys;
end
