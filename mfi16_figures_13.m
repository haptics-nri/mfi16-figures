% part 13 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
           %%
            fprintf('Loading data for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            dataset = [materials{mi} reps{ri} tools{ti}];
            prefix = [DATADIR filesep date{ri} filesep dataset filesep];
            
            [v{mi,ri,ti}, int{mi,ri,ti}, acc{mi,ri,ti}] = load_stick(prefix);
