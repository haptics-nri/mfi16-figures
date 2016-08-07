% part 38 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
                confusion{cvi} = zeros(5,5);
                predictions{cvi} = svmpredict(zeros(size(val_vectors,1),1), val_vectors(:,2:end), model, '-q');

                for j=1:5
                    choose_stick_39
                end

                cv_acc(cvi) = sum(diag(confusion{cvi}))/sum(sum(confusion{cvi}));
                fprintf('\tFold %d: MC %g%%\n', cvi, 100*cv_acc(cvi));
