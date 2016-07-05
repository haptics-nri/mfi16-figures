% part 17 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    fprintf('%d calculating pre-features... ', i);
    j = endeff_idx(eps(i).endeff);
    if isfield(eps(i).data, 'f_comp')
        choose_stick_18
    end
    fprintf('done\n');
