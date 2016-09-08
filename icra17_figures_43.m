% part 43 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    figure(fv1);
    subplot(5,1,i);
    idx = cell2mat(train_features(:,1))==i;
    g = f(idx,:);
    imagesc(g);%, [allmin allmax]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box off; axis off;
    text(0.3, size(g,1)/2, ...
         material_names{i}, ...
         'FontSize', 14, ...
         'HorizontalAlignment', 'right', ...
         'Interpreter', 'tex');
     
    %figure(fv2);
    %subplot(1,6,i);
    %cor = corrcoef(g(:,1:end-6));
    %imagesc(cor);
