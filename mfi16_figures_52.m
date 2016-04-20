% part 52 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% feature vectors -- first set gsi to optimal and run the grid search iter
fv1 = figure;
fv2 = figure;
f = romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode);
m = mean(f);
for i=1:5
    mfi16_figures_53
end
figure(fv1);
colormap jet;
axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
set(colorbar('ticks',[]), 'edgecolor','none');
print -dpdf mfi16_feature_vectors.pdf;
figure(fv2);
subplot(1,6,6);
colormap jet;
colorbar;

