% part 13 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
        % replace NUC timestamps with Teensy deltas
        ep = data(materials{m});
        ep.nuc_t = ep.int(:,1);
        ep.int(:,1) = ep.int(1,1) + cumsum(bitand(ep.dt, 65535))/1e6;
        data(materials{m}) = ep;
