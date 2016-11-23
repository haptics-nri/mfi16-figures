% part 19 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
    if dosub
        icra17_figures_20
    else
        icra17_figures_21
    end
    set(gca, 'FontSize',12);
    hold on;
    for j=0:150:1500
        icra17_figures_22
    end
    datums = plot(d{i}.iws(d{i}.a:d{i}.b,1)-d{i}.iws(1,1), ...
         [...
          d{i}.iws(d{i}.a:d{i}.b,4) ...
          sqrt(sum(d{i}.iws(d{i}.a:d{i}.b,2:3).^2,2)) ...
          dft321(d{i}.ai(d{i}.a:d{i}.b,2:4))*10 ...
          sqrt(sum(d{i}.vel(d{i}.a:d{i}.b,:).^2,2))/10 ...
         ]);
    axis([d{i}.iws(d{i}.a,1)-d{i}.iws(1,1) d{i}.iws(d{i}.b,1)-d{i}.iws(1,1) -15 30]);
    set(gca, 'XTick', round(d{i}.iws(d{i}.a + (225:450:1499)',1) - d{i}.iws(1,1), 2));
    legend(datums([1 2 4 3]), 'Normal force (N)', 'Tangential force (N)', 'Tip speed (cm/s)', 'Acceleration (cm/s^2)', 'location','northeast');
    xlabel('Time (s)');
    box on;
    if ~dosub
        icra17_figures_24
    end
    
    if dosub
        icra17_figures_25
    else
        icra17_figures_26
    end
    
    posted = d{i}.posted;
    posted = posted';
    
    cla;
    imagesc(posted, clim);
    set(gca, 'FontSize',12);
    colormap jet;
    box off;
    sub.XTick = 1:10;
    sub.YTickLabel = [];
    xlabel('Feature vectors');
    sub.XRuler.Axle.Visible = 'off';
    sub.YRuler.Axle.Visible = 'off';
    labels = repmat({''}, [1 5]);
    for j=1:5
        icra17_figures_27
    end
    things = {'Fn', 'Ft', 'V'};
    for thing=1:length(things)
        icra17_figures_28
    end
    for j=1:length(labels)
        icra17_figures_29
    end
    if ~dosub
        icra17_figures_30
    end
