% part 43 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        subplot(5,1,j);
        g = f(cell2mat(train_features(:,1))==j, :);
        g = bsxfun(@minus, g, m);
        g = bsxfun(@rdivide, g, max(g) - min(g));
        g = [g mean(g(:,[end-5 end-3 end-1]), 2)];
        g = sortrows(g, size(g,2));
        g = g(:,1:end-1);
        imagesc(g);
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        box off; axis off;
        text(0.3, size(g,1)/2, ...
             material_names{j}, ...
             'FontSize', 14, ...
             'HorizontalAlignment', 'right', ...
             'Interpreter', 'tex');
