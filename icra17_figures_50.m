% part 50 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
        if i == j
            icra17_figures_51
        else
            icra17_figures_52
        end
        text(j, i, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(:,j))), ...
             'FontSize',14, 'Color',c, 'horizontalalignment','center');
