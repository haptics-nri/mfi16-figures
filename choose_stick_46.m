% part 46 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
            endeffs(i).confusion(j,k) = nnz(endeffs(i).predictions(test_vectors(:,1)==j) == k);
