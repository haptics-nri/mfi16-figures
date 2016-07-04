% part 13 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    figs(e) = figure;
    for s=1:length(order)
        choose_stick_14
    end
    suplabel(eps(k).endeff, 't');
    print('-depsc', sprintf('choose_stick__%d_eighths.eps', eval(eps(k).endeff)*8));
    order = fliplr(order);
