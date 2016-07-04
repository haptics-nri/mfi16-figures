% part 36 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    endeffs(i).model = svmtrain(endeffs(i).vectors.train(:,1), endeffs(i).vectors.train(:,2:end), '-q -m 1000 -s 1 -t 2 -n 0.1 -g 7.0');
    endeffs(i).predictions = svmpredict(zeros(size(endeffs(i).vectors.test,1),1), endeffs(i).vectors.test(:,2:end), endeffs(i).model, '-q');
    endeffs(i).confusion = zeros(5,5);
    for j=1:5
        choose_stick_37
    end
    endeffs(i).accuracy = sum(diag(endeffs(i).confusion))/sum(sum(endeffs(i).confusion));
