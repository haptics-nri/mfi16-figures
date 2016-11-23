% part 31 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    figure;
    set(colorbar, 'XTick', [0 1], 'XTickLabel', {'min' 'max'});
    colormap jet;
    print -dpdf icra17_fv_colorbar.pdf;
