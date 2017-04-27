function errors = lint(episode)

    checks = {
        'sensor rates'      @sensor_rates
        'sensor magnitudes' @sensor_mags
        'motion track'      @motion_track
        'images'            @images
    };

    errors = {};
    for i = 1:size(checks,1)
        fprintf('Checking %s...\n', checks{i,1});
        errors = [errors
                  checks{i,2}(episode)];
    end

end

function check(msg, val, err)
    fprintf('\tChecking %s...\n', msg);
    if ~val
        assignin('caller', 'errors', [evalin('caller', 'errors'); err]);
    end
end

function errors = sensor_rates(episode)
    rates = {
        'v'        15
        'int'      3000
        'acc'      3000
        'imu.acc'  1340
        'imu.gyro' 750
        'imu.mic'  3000

        'bvei'     3000
        'biws'     3000
        'bai'      3000
    };
    if strcmp(episode.flow, 'optocam')
        rates(end+1,:) = {'sensor.opto'    1000};
    elseif strcmp(episode.flow, 'biocam')
        rates(end+1,:) = {'sensor.bio'     100 };
    end

    errors = {};

    for i = 1:size(rates,1)
        k = rates{i,1};

        if strfind(k, '.')
            kk = strsplit(k, '.');
            rate = 1/mean(diff(episode.(kk{1}).(kk{2})(:,1)));
        else
            rate = 1/mean(diff(episode.(k)(:,1)));
        end

        check(sprintf('sensor rate of %s', k), ...
              abs(rates{i,2} - rate)/rates{i,2} <= 0.05, ...
              sprintf('Sensor rate of %s is %g but it should be %g!', k, rate, rates{i,2}));
    end
end

function errors = sensor_mags(episode)
    errors = {};

    check('force', ...
          all(episode.biws(:,4) > 0), ...
          'Rectified Z force should be all positive!');

    if strcmp(episode.flow, 'stickcam')
        mic_freq = abs(fft(episode.imu.mic(:,2) - mean(episode.imu.mic(:,2))));
        check('sound', ...
              max(mic_freq(10:round(end/4))) / max(mic_freq(round(end/4):round(end/2))) < 2, ...
              'Microphone frequency content is too uneven!');
    end

    dacc = bsxfun(@minus, episode.bai(:,2:4), mean(episode.bai(:,2:4)));
    check('digital acceleration magnitude', ...
          all(std(dacc) > 0.1), ...
          'End-effector accelerometer data should be nonzero!');

    aacc = bsxfun(@minus, episode.imu.acc(:,2:4), mean(episode.imu.acc(:,2:4)));
    check('analog acceleration', ...
          all(std(aacc) > 0.1), ...
          'Analog accelerometer data should be nonzero!');

    agyro = bsxfun(@minus, episode.imu.gyro(:,2:4), mean(episode.imu.gyro(:,2:4)));
    check('analog gyroscope', ...
          all(std(agyro) > 0.1), ...
          'Analog gyroscope data should be nonzero!');
end

function errors = motion_track(episode)
    errors = {};

    % FIXME why are there any NaNs at all
    first_nan = min([...
        find(isnan(episode.bvei(:,2)),1) ...
        find(isnan(episode.bvei(:,3)),1) ...
        find(isnan(episode.bvei(:,4)),1) ...
        ]);
    
    check('position NaN', ...
          first_nan/size(episode.bvei, 1) > 0.98, ...
          'Position should not be NaN until the very end!');

    p = RANSAC(episode.bvei(1:first_nan-1,2:4)', ...
               struct('sigma', 0.5, ...
                      'est_fun', @estimate_plane, ...
                      'man_fun', @error_plane, ...
                      'verbose', false));
    normal = p.Theta(1:3)/sqrt(sum(p.Theta(1:3).^2));

    check('position in plane', ...
          nnz(p.CS)/first_nan > 0.5, ...
          'Motion should be in-plane!');

    check('plane orientation', ...
          abs(dot(normal, [0 0 1])) > 0.98, ...
          'Motion track should in the horizontal plane!');

    plane = fitPlane(episode.bvei(p.CS,2:4));
    inplane = planePosition(episode.bvei(p.CS,2:4), plane);
    check('exploration area', ...
          nth(2, 2, @convhull, inplane) > 50^2, ...
          'Exploration area should be at least 50 mm^2!');
end

function errors = images(episode)
    errors = {};

    frames = episode.april.keys;
    img = imread(fullfile(episode.datadir, episode.date, episode.flow, episode.i, 'bluefox', sprintf('bluefox%d.png', frames{1})));
    
    check('image brightness', ...
          mean(img(:)) > 80, ...
          'Image is too dark!');
end

