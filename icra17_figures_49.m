% part 49 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% confusion matrices -- first set gsi to optimal and run the test set
figure;
imagesc(bsxfun(@rdivide, mc_test_confusion, sum(mc_test_confusion, 2)), [0 1]);
colormap(flipud(gray));
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [1 2 3 4 5];
ax.XTickLabel = {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'};
ax.YTickLabel = {'ABS', 'glitter paper', 'silk', 'vinyl', 'wood'};
ax.FontSize = 14;
xlabel('Detected material');
ylabel('Actual material');
ax.XLabel.Position = ax.XLabel.Position + [0 0.1 0];
for i=1:length(materials)
    icra17_figures_50
end
print -dpdf icra17_confusion_precision.pdf;

fprintf('\n');
fprintf('Surface & Precision & Recall & $F_1$ score\n');
for i=1:length(materials)
    icra17_figures_54
end

