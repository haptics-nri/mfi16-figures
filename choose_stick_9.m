% part 9 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    fprintf('%d\n', i);
    if std(diff(eps(i).data.force(100:end,1))) < 2e-4
        choose_stick_10
    end
