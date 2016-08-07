% part 48 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    %[~,gsi] = max(endeffs(i).gs_acc);
    gs_nbins = 7;%3;%nbins(gs_idx(gsi,1));
    gs_binmode = 'perceptual';%binmode{gs_idx(gsi,2)};
    gs_alpha = 0.2;%alpha(gs_idx(gsi,3));
    gs_nu = 0.1;%nu(gs_idx(gsi,4));
    gs_gamma = 7;%gamma(gs_idx(gsi,5));
    gs_stmode = true;%stmode(gs_idx(gsi,6));
    
    figure;
    f = romano_features('post', endeffs(i).features(:,2:end), gs_nbins, gs_binmode, gs_alpha, gs_stmode);
    m = mean(f);
    for j=1:5
        choose_stick_49
    end
    
    a = subplot(5,1,5);
    axis on;
    a.XRuler.Axle.Visible = 'off';
    a.YRuler.Axle.Visible = 'off';
    labels = {};
    for j=1:gs_nbins
        choose_stick_50
    end
    things = {'F_N', 'V', 'F_T'};
    for thing=1:length(things)
        choose_stick_51
    end
    for j=1:length(labels)
        choose_stick_53
    end
    colormap jet;
    axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
    set(colorbar('ticks',[]), 'edgecolor','none');
    text(1.08, 1, '1', 'fontsize',14);
    text(1.08, 0, '0', 'fontsize',14);
    
    subplot(5,1,1);
    title(sprintf('endeff %d feature vectors', i));
