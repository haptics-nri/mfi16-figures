% part 38 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
            endeffs(i).confusion(j,k) = nnz(endeffs(i).predictions(endeffs(i).vectors.test(:,1)==j) == k);
