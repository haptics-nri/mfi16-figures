% part 51 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
        if i == j
            icra17_figures_52
        else
            icra17_figures_53
        end
        text(j, i, sprintf('%.3f', mc_test_confusion(i,j)/sum(mc_test_confusion(i,:))), ...
             'FontSize',14, 'Color',c, 'HorizontalAlignment','center');