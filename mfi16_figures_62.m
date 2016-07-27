% part 62 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
        if i == j
            mfi16_figures_63
        else
            mfi16_figures_64
        end
        text(j, i, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(:,j))), ...
             'FontSize',14, 'Color',c, 'horizontalalignment','center');
