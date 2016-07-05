% part 20 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    j = endeff_idx(eps(i).endeff);
    if isfield(eps(i).data, 'f_comp') && ~strcmp(eps(i).material, 'weave')
        choose_stick_21
    end
