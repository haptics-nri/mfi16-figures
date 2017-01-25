% Plots a confusion matrix in a nice-looking figure
% Clobbers the current figure.
% inputs:
%   - confusion: the NxN confusion matrix
%   - materials: 1D cell array with class label strings
%   - fontsize, fontface, fontrot: font size, name and rotation of the class labels (axis labels are always 16pt Arial)
%   - textpct: boolean flag, whether to print the confusion values on the grid
function fig_confusion(confusion, materials, fontsize, fontface, fontrot, textpct)
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
end

