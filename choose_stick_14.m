% part 14 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        subplot(2,3,s);
        k = (e-1)*30 + (order(s)-1)*5 + 1;
        if ~isfield(eps(k).data, 'f_comp')
            choose_stick_15
        end
        a = eps(k).peaks(2) + 1500;
        b = eps(k).peaks(3) - 1500;
        plot(eps(k).data.f_comp(a:b,1)-eps(k).data.f_comp(1,1), eps(k).data.f_comp(a:b,2:4), ...
             eps(k).data.v_end(a:b,1)-eps(k).data.f_comp(1,1), bsxfun(@minus, eps(k).data.v_end(a:b,2:3), eps(k).data.v_end(a,2:3))/10 + 20);
        title(eps(k).material);
        axis tight;
