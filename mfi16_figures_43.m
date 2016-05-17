% part 43 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
            % evaluate by comparing all OCSVMs and the MCSVM
            %oc_confusion{cvi} = zeros(length(materials));
            mc_confusion{cvi} = zeros(length(materials));
            %prob = zeros(size(val_vectors,1),length(materials));
            %for mi=1:length(materials)
            %    prob(:,mi) = rabaoui_dissim(models{mi}, val_vectors(:,2:end));
            %end
            %[~, oc_answers{cvi}] = min(prob, [], 2);
            mc_answers{cvi} = svmpredict(zeros(size(val_vectors,1),1), val_vectors(:,2:end), models{end}, '-q');

            for i=1:length(materials)
                mfi16_figures_44
            end

            cv_acc(cvi) = sum(diag(mc_confusion{cvi}))/sum(sum(mc_confusion{cvi}));
            fprintf('\tFold %d: MC %g%%\n', cvi, 100*cv_acc(cvi));
