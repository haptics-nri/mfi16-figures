% part 54 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% confusion matrices -- first set gsi to optimal and run the test set
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 1)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.YTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.FontSize = 14;
xlabel('Detected material');
ylabel('Actual material');
ax.XLabel.Position = ax.XLabel.Position + [0 0.1 0];
for i=1:length(materials)
    mfi16_figures_55
end
print -dpdf mfi16_confusion_precision.pdf;
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 2)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.YTickLabel = {'ABS', 'paper plate', 'folder', 'MDF', 'canvas'};
ax.FontSize = 14;
xlabel('Detected material');
ylabel('Actual material');
ax.XLabel.Position = ax.XLabel.Position + [0 0.1 0];
for i=1:length(materials)
    mfi16_figures_59
end
print -dpdf mfi16_confusion_recall.pdf;
