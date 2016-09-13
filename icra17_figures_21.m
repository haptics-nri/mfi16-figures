% part 21 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
   %%
    features = [features
                num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
            
    new_feats = romano_features('pre', ep.biws, ep.bvei, ep.bai, mass, 150, [5 .5], ep.bss);
    bfeatures = [bfeatures
                num2cell(repmat(m, size(new_feats,1), 1)) new_feats];
