% Test a set of hyperparameters for WHC 17
function [confusions, accs, models] = whc17_test(materials, win, features, labels, cv, mode, grid_values, grid_params, gs_acc, filename)
    confusions = zeros(cv.NumTestSets, length(materials), length(materials));
    accs = zeros(1, cv.NumTestSets);
    models = cell(1, length(materials));

    for cvi=1:cv.NumTestSets
        [~,gsi] = max(gs_acc(cvi,:));
        gp = grid_params(gsi,:);
        gs_nbins   = grid_values{1}(gp(1));
        gs_binmode = grid_values{2}{gp(2)};
        gs_alpha   = grid_values{3}(gp(3));
        gs_nu      = grid_values{4}(gp(4));
        gs_gamma   = grid_values{5}(gp(5));

        fprintf('%s K=%d win=%d: n=%d bins=%s α=%g ν=%g ɣ=%g\n', mode, cvi, win, gs_nbins, gs_binmode, gs_alpha, gs_nu, gs_gamma);
        switch mode
            case 'romano'
                tr_rf = romano_features('post', features{1}(cv.training(cvi),:), gs_nbins, gs_binmode, gs_alpha, 1);
                te_rf = romano_features('post', features{1}(cv.test    (cvi),:), gs_nbins, gs_binmode, gs_alpha, 1);
            case 'steinbach'
                tr_sf = features{2}(cv.training(cvi),:);
                te_sf = features{2}(cv.test    (cvi),:);
            case 'both'
                tr_rf = romano_features('post', features{1}(cv.training(cvi),:), gs_nbins, gs_binmode, gs_alpha, 1);
                te_rf = romano_features('post', features{1}(cv.test    (cvi),:), gs_nbins, gs_binmode, gs_alpha, 1);
                tr_sf = features{2}(cv.training(cvi),:);
                te_sf = features{2}(cv.test    (cvi),:);
        end
        switch mode
            case 'romano'
                train_features = tr_rf;
                test_features  = te_rf;
            case 'steinbach'
                train_features = tr_sf;
                test_features  = te_sf;
            case 'both'
                train_features = [tr_rf tr_sf];
                test_features  = [te_rf te_sf];
        end
        train_labels = labels(cv.training(cvi),:);
        test_labels  = labels(cv.test    (cvi),:);
        
        fprintf('Testing...\n');

        [confusions(cvi,:,:), accs(cvi), models{cvi}] = mc_svm(materials, train_labels, train_features, test_labels, test_features, gs_nu, gs_gamma);
        fprintf('Test set accuracy (K=%d, win=%d, r): %g\n', cvi, win, accs(cvi));
    end

    fprintf('Plotting...\n');
    conf = figure;
    fig_confusion(squeeze(mean(confusions)), materials, 10, 'Courier', 90, 0, false);
    grid on;
    print('-dpdf', filename);
end

