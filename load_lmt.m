% loads data from the LMT DB
% 
% parameters:
% - datadir: directory containing data
% - traintest: 'train' to load an episode from the training set, 'test' for testing set
% - label: 'surface label, e.g. 'G6Carpet' or 'G9TextileVersion1'
% - id: trial number or name to select the episode, e.g. '0' (if traintest == 'train') or 'Andrea' (if traintest == 'test')
% returns:
% - accel: contact vibration, Nx1 (combined by DFT321)
% - friction: tangential force, Nx1
% - sound: sound amplitude, Nx2
function [accel, friction, sound] = load_lmt(datadir, traintest, label, id)

    if strcmp(traintest, 'train')
        accel = load_1d(fullfile(datadir, 'Training', 'Accel'), label, sprintf('query%s', id));
        friction = load_1d(fullfile(datadir, 'Training', 'Friction'), label, 'continuosFric');
        sound = load_wav(fullfile(datadir, 'Training', 'Sound'), label, sprintf('sound%s', id));
    else
        accel = load_1d(fullfile(datadir, 'Testing', 'AccelScansDFT321'), label, id);
        friction = load_1d(fullfile(datadir, 'Testing', 'FricScans'), sprintf('%s_Friction', label), id);
        sound = load_wav(fullfile(datadir, 'Testing', 'SoundScans'), sprintf('%s_Sound', label), id);
    end

end

function signal = load_1d(path, label, id)
    file = fullfile(path, sprintf('%s_%s.txt', label, id));
    %fprintf('\tloading %s...\n', file);
    data = csvread(file);
    signal = struct('data', data, 'Fs', 10000);
end

function signal = load_wav(path, label, id)
    file = fullfile(path, sprintf('%s_%s.wav', label, id));
    %fprintf('\tloading %s...\n', file);
    [snd, fs] = audioread(file);
    signal = struct('data', snd, 'Fs', fs);
end

