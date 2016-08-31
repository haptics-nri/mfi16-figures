% part 13 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    ep = data(materials{m});
    
    [~,~,~,~,~,~, vei, ai, ~,~,~, iws] = process_stick(ep.v, ep.int, ep.acc, mass, [0;0;0], H_vic2bod, H_m402bod, H_bal2imu, -ep.off);
    ep.vei = vei;
    ep.ai = ai;
    ep.iws = iws;
    
    data(materials{m}) = ep;
