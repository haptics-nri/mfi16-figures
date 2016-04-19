% part 26 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
           %%
            fprintf('Extracting features for %s on %s material, rep #%s\n', tools{ti}, materials{mi}, reps{ri});
            
            fprintf('\tfriction coefficient\n');
            [mu, err] = extract_friction(intworldsub{mi,ri,ti}, vendint{mi,ri,ti});
            mu_k{mi,ri,ti} = real([mu err]);
            
            fprintf('\tspringiness\n');
            [k, z0, err] = extract_springiness(intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, mass);
            spring{mi,ri,ti} = [k z0 err];
            
            fprintf('\tvibration power\n');
            power{mi,ri,ti} = extract_vibration(intworldsub{mi,ri,ti}, vendint{mi,ri,ti}, accint{mi,ri,ti}, mass);
