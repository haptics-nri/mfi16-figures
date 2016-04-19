% part 31 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
           %%
            new_feats = romano_features('pre', intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, accint{mi,ri,ti}, mass, 0.05, [20 3]);
                                                                         % FIXME reexamine these thresholds, use forcefiltsub(:,3) instead of forcemag
