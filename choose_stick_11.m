% part 11 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
            j = endeff_idx(eps(i).endeff);
            [~,~,~,~,~,~, eps(i).data.v_end, eps(i).data.a_end, ~,~,~, eps(i).data.f_comp] = process_stick(eps(i).data.vicon, eps(i).data.force, eps(i).data.acc, masses(j), coms(j,:)', H_vic2bod, H_m402bod, H_bal2imu, -eps(i).offset);
