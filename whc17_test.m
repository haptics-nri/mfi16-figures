% Test a set of hyperparameters for WHC 17
function [confusion, model] = whc17_test(materials, wins, features, labels, split_idx, w, mode, nbins, binmode, alpha, nu, gamma, filename)
    fprintf('Calculating features...\n');
    switch mode
        case 'romano'
            tr_rf = romano_features('post', features{w}{1}(split_idx{w}==1,:), nbins, binmode, alpha, 1);
            te_rf = romano_features('post', features{w}{1}(split_idx{w}==2,:), nbins, binmode, alpha, 1);
        case 'steinbach'
            tr_sf = features{w}{2}(split_idx{w}==1,:);
            te_sf = features{w}{2}(split_idx{w}==2,:);
        case 'both'
            tr_rf = romano_features('post', features{w}{1}(split_idx{w}==1,:), nbins, binmode, alpha, 1);
            te_rf = romano_features('post', features{w}{1}(split_idx{w}==2,:), nbins, binmode, alpha, 1);
            tr_sf = features{w}{2}(split_idx{w}==1,:);
            te_sf = features{w}{2}(split_idx{w}==2,:);
    end
    switch mode
        case 'romano'
            train_features = tr_rf;
            test_features = te_rf;
        case 'steinbach'
            train_features = tr_sf;
            test_features = te_sf;
        case 'both'
            train_features = [tr_rf tr_sf];
            test_features = [te_rf te_sf];
    end
    train_labels = labels{w}(split_idx{w}==1,:);
    test_labels = labels{w}(split_idx{w}==2,:);
    
    fprintf('Testing...\n');

    [confusion, acc, model] = mc_svm(materials, train_labels, train_features, test_labels, test_features, nu, gamma);
    fprintf('Test set accuracy (win=%d, r): %g\n', wins(w), acc);

    fprintf('Plotting...\n');
    conf = figure;
    fig_confusion(confusion, materials, 10, 'Courier', 90, false);
    grid on;
    print('-dpdf', filename);
end

