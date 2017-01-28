% Plots a confusion matrix in a nice-looking figure
% Clobbers the current figure.
% inputs:
%   - confusion: the NxN confusion matrix
%   - materials: 1D cell array with class label strings
%   - fontsize, fontface, fontrot: font size, name and rotation of the class labels (axis labels are always 16pt Arial)
%   - xloff: amount to kick down the xlabel
%   - textpct: boolean flag, whether to print the confusion values on the grid
function fig_confusion(confusion, materials, fontsize, fontface, fontrot, xloff, textpct)
    clf;
    imagesc(bsxfun(@rdivide, confusion, sum(confusion, 2)), [0 1]);
    axis square;
    colormap(flipud(gray));
    ax = gca;
    ax.XTick = 1:length(materials);
    ax.XTickLabel = materials;
    ax.XTickLabelRotation = fontrot;
    ax.YTick = 1:length(materials);
    ax.YTickLabel = materials;
    ax.FontSize = fontsize;
    ax.XLabel.Position(2) = ax.XLabel.Position(2) + xloff;

    if textpct
        for i=1:length(materials)
            for j=1:length(materials)
                if i == j
                    c = 'white';
                else
                    c = 'black';
                end
                text(j, i, sprintf('%.3f', confusion(i,j)/sum(confusion(i,:))), ...
                     'FontSize',14, 'Color',c, 'HorizontalAlignment','center');
            end
        end
    end

    xlabel('Detected material', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Actual material', 'FontName', 'Arial', 'FontSize', 16);
    set(gca, 'FontName', fontface);
    axis normal;
    
    fprintf('\n');
    fprintf('Surface & Precision & Recall & $F_1$ score \\\\\n\\hline\n');
    prec = zeros(1,length(materials));
    rec = zeros(1,length(materials));
    f1 = zeros(1,length(materials));
    for i=1:length(materials)
        fprintf('%s  &  ', materials{i});
        others = setdiff(1:5, i);
        tp = confusion(i,i);
        fp = sum(confusion(others,i));
        fn = sum(confusion(i,others));
        prec(i) = tp/(tp+fp);
        rec(i) = tp/(tp+fn);
        f1(i) = 2*prec(i)*rec(i)/(prec(i)+rec(i));
        fprintf('%.3f  &  %.3f  &  %.3f', prec(i), rec(i), f1(i));
        fprintf(' \\\\ \n');
    end
    fprintf('\\emph{mean}  &  %.3f  &  %.3f  &  %.3f\n', mean(prec), mean(rec), mean(f1));
end

