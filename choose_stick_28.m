% part 28 of /Users/alex/Documents/research/proton/code/calibration/motion/choose_stick.m
    % CV partition
    cv = cvpartition(cell2mat(endeffs(i).features(endeffs(i).split==1, 1)), 'KFold', 5);
    confusion = cell(1, cv.NumTestSets);
    predictions = cell(1, cv.NumTestSets);
    endeffs(i).gs.partition = cv;
    endeffs(i).gs.cv_acc = cell(size(gs_idx,1),1);
    endeffs(i).gs.gs_acc = zeros(size(gs_idx,1),1);
    
    gs_acc = zeros(size(gs_idx,1),1);
    clear romano_features; % clear persistent vars
    
    elapsed = tic;
    for gsi=1:size(gs_idx,1)
        choose_stick_29
    end
    endeffs(i).gs_acc = gs_acc;
