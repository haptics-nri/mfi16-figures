% part 54 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% feature vectors -- first set gsi to optimal and run the grid search iter
fv1 = figure;
%fv2 = figure;
f = romano_features('post', train_features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode);
m = mean(f);
for i=1:5
    mfi16_figures_55
end
figure(fv1);
a = subplot(5,1,5);
axis on;
a.XRuler.Axle.Visible = 'off';
a.YRuler.Axle.Visible = 'off';
%a.XTick = 1:10;
labels = {};
for i=1:gs_nbins
    mfi16_figures_56
end
things = {'F_N', 'V', 'F_T'};
for thing=1:length(things)
    mfi16_figures_57
end
%a.XTickLabels = labels;
%a.XTickLabelRotation = 70;
%a.TickLength = [0 0];
for i=1:length(labels)
    mfi16_figures_59
end
colormap jet;
axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
set(colorbar('ticks',[]), 'edgecolor','none');
text(1.08, 1, '1', 'fontsize',14);
text(1.08, 0, '0', 'fontsize',14);
print -dpdf mfi16_feature_vectors.pdf;
%figure(fv2);
%subplot(1,6,6);
%colormap jet;
%colorbar;

