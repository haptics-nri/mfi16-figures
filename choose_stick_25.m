% part 25 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
        confusion = zeros(5,5);
        predictions = svmpredict(zeros(size(test_vectors,1),1), test_vectors(:,2:end), model, '-q');

        for j=1:5
            choose_stick_26
        end

        fprintf('\tEndeff %d: MC %g%%\n', i, sum(diag(confusion))/sum(sum(confusion)));
