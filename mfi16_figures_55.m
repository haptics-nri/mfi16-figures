% part 55 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
    figure(fv1);
    subplot(5,1,i);
    g = f(cell2mat(train_features(:,1))==i, :);
    g = bsxfun(@minus, g, m);
    g = bsxfun(@rdivide, g, max(g) - min(g));
    g = [g mean(g(:,(end-6):end), 2)];
    g = sortrows(g, size(g,2));
    g = g(:,1:end-1);
    imagesc(g);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box off; axis off;
    text(0.3, size(g,1)/2, ...
         material_names{i}, ...
         'FontSize', 14, ...
         'HorizontalAlignment', 'right', ...
         'Interpreter', 'tex');
     
    figure(fv2);
    subplot(1,6,i);
    cor = corrcoef(g(:,1:end-6));
    imagesc(cor);
