function gs_plots(acc, idx, vars)
    figure;
    for i=1:size(vars,1)
        subplot(size(vars,1), 1, i);
        plot(vars{i,2}(idx(:,i)), acc, '.');
        hold on;
        plot(vars{i,2}, grpstats(acc, idx(:,i), @median), 'r.-', 'MarkerSize',30);
        hold off;
        xlabel(vars{i,1});
        ylabel('Accuracy');
        title(sprintf('Accuracy vs %s', vars{i,1}));
    end
end
