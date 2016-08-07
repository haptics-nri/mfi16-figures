% part 18 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        % FIXME change this back to use analog accelerometers
        %eps(i).features.pre = romano_features('pre', eps(i).data.f_comp, eps(i).data.v_end, eps(i).data.a_end, masses(j), 0.05, [20 3], [eps(i).peaks(2)+1000 eps(i).peaks(3)-1000]);
        % use high-passed force instead of accelerometer
        filtered = filter([1 .02-1], [.02-1 0], eps(i).data.f_comp(:,2:4));
        eps(i).features.pre = romano_features('pre', eps(i).data.f_comp, eps(i).data.v_end, filtered, masses(j), 0.05, [20 3], [eps(i).peaks(2)+1000 eps(i).peaks(3)-1000]);
        eps(i).features.post = romano_features('post', eps(i).features.pre, 3, 'perceptual', 0.2, true);
