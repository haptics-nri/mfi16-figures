function [features, split_idx, bfeatures, bsplit_idx] = icra17_svm(data, mass)

    materials = data.keys;

    features = cell(0, 5); % cols: label, vibration, speed, normal, tangential
    bfeatures = cell(0, 5); % cols: label, vibration, speed, normal, tangential
    for m = 1:length(materials)
        ep = data(materials{m});

        fprintf('Romano features for %s\n', materials{m});
        %%
        new_feats = romano_features('pre', ep.iws, ep.vei, ep.ai, mass, 150, [5 .5], ep.ss);
        %%
        features = [features
                    num2cell(repmat(m, size(new_feats,1), 1)) new_feats];

        new_feats = romano_features('pre', ep.biws, ep.bvei, ep.bai, mass, 150, [5 .5], ep.bss);
        bfeatures = [bfeatures
                    num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
    end

    % test/train split

    % 4/5 train, 1/5 test
    split_idx = randsample(1:2, size(features,1), true, [4/5 1/5]);
    bsplit_idx = randsample(1:2, size(bfeatures,1), true, [4/5 1/5]);

end
