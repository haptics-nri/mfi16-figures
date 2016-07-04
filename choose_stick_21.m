% part 21 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        endeffs(j).features = [endeffs(j).features
                               num2cell(repmat(material_idx(eps(i).material), [size(eps(i).features.pre,1) 1])) eps(i).features.pre];
