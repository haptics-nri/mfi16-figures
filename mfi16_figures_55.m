% part 55 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
        if i == j
            mfi16_figures_56
        end
        text(i, j, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(i,:))), ...
             'Color',c, 'HorizontalAlignment','center');
