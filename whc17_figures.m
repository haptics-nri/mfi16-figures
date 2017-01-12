%% generates figures for the WHC 2017 paper

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

DATADIR = '/Volumes/shared/Projects/Proton Pack/Data';

%% sphere calibration (see go_sphere_again.m)

% vicon data
x = csvload([DATADIR '/20160726/sphere/4/vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

[c, r, e, cs] = sphereFit_ransac(x(:,2:4), 50); % FIXME allowing a huge amount of noise to get a reasonable number of inliers

d = nan([size(x,1) 3]);
for i = find(cs)
    cc = xfconv(x(i,2:7)) \ [c 1]';
    d(i,:) = cc(1:3);
end
d = nanmean(d);

% the product of sphere calibration is d
spherecalib.x = x;
spherecalib.c = c;
spherecalib.r = r;
spherecalib.cs = cs;
spherecalib.d = d;

%% free-space calibration (see go_bias.m)

date = '20160226';
offset = 22.848; % manual (from taps)
material = 'free';
tool = 'stick';
rep = '1';

set = [material rep tool];
prefix = [DATADIR filesep date filesep set filesep];

v = csvload([prefix 'vicon.tsv'], ...
            {'Timestamp', ...
             'proton_Root_T_X_', 'proton_Root_T_Y_', 'proton_Root_T_Z_', ...
             'proton_Root_A_X_', 'proton_Root_A_Y_', 'proton_Root_A_Z_'}, ...
            {'Delimiter', '\t'});

int = csvload([prefix 'teensy.ft.csv'], ...
              [{'Timestamp'}, ...
               arrayfun(@(x) ['FT' num2str(x)], 0:29, 'UniformOutput',false)]);
           
int = process_mini40(100*ones(size(int,1),1), int, zeros(1,6));

a = round(size(int,1)*1/5);
b = round(size(int,1)*4/5);
int_s = int(a:b,:);

[~,aa] = min(abs(v(:,1)-int(a,1)));
[~,bb] = min(abs(v(:,1)-int(b,1)));
v_s = v(aa:bb,:);

[c,r,e,inliers] = sphereFit_ransac(int_s(:,2:4)); % 94% inliers

[~, start] = min(abs(int_s(:,1) - v(1,1)));

R_m402bod = [ 0       0       1
              1       0       0
              0       1       0 ];

grav = zeros(size(int_s,1), 4);
grav(:,1) = int_s(:,1) - offset;
ideal_grav = grav;
for tm=1:size(grav,1)
    [delta, tv] = min(abs(v(:,1) - grav(tm,1)));
    assert(delta < 0.02);
    Rvw = xfconv(v(tv,5:7));
    grav(tm,2:4) = R_m402bod * int_s(tm,2:4)';
    ideal_grav(tm,2:4) = Rvw\[0 0 -r]';
end

[R, t, err] = rigid_ransac(ideal_grav(inliers,2:4), grav(inliers,2:4), 10); % 100% inliers??

% the product of free-space calibration is R
freecalib.int_s = int_s;
freecalib.grav = grav;
freecalib.ideal_grav = ideal_grav;
freecalib.c = c;
freecalib.r = r;
freecalib.cs = inliers;
freecalib.R = R;
freecalib.t = t;

%% end-effector weighing

date = '20160811';
[~,f] = load_stick([DATADIR filesep date filesep 'weigh/1/']);
[c,r] = sphereFit_ransac(f(:,2:4)); % 98.9% inliers

mass = r/9.81;

%% process calibrations: use d+R to calculate x+y+z for H_vic2bod + H_bal2imu

H_vic2bod = [ 0.9912 -0.0238 -0.1302    0
              0.0164  0.9982 -0.0575   36
              0.1314  0.0549  0.9898 -511.32
              0       0       0         1   ];
H_m402bod = [ 0       0       1       108.99
              1       0       0         0.53
              0       1       0        -2.98
              0       0       0         1   ];
H_bal2imu = [ 1       0       0       254.84
              0       1       0         0
              0       0       1         0
              0       0       0         1   ];

syms x y z real;
R = freecalib.R';
H1 = [R [0 y z]';      0 0 0 1];
H2 = [eye(3) [x 0 0]'; 0 0 0 1];
S = solve(H1 * H2 * [0 0 0 1]' == [d 1]');

H_vic2bod(1:3,1:3) = R;
H_vic2bod(2,4) = S.y;
H_vic2bod(3,4) = S.z;
H_bal2imu(1,4) = S.x;

%% load data

% dataset parameters
date = '20170109'; % FIXME also 20161229
flowtype = 'stickcam';

[data1, materials1] = icra17_load(DATADIR, '20161229', flowtype, @(x) false);
[data2, materials2] = icra17_load(DATADIR, '20170109', flowtype, @(x) false);
materials = [materials1 materials2];
data = containers.Map;
for m=1:length(materials1)
    data(materials1{m}) = data1(materials1{m});
end
for m=1:length(materials2)
    data(materials2{m}) = data2(materials2{m});
end
clear data1 data2
data = icra17_process('bluefox', data, mass, H_vic2bod, H_m402bod, H_bal2imu);

%% calculate features

for m=1:length(materials)
    fprintf('Calculating %d/%d %s... ', m, length(materials), materials{m});
    
    d = data(materials{m});
    
    % how to find a tap
    % 1. take median of Z position
    % 2. find first point where Z position goes 5mm above median
    % 3. tap starts at the first point after that where the Z derivative goes negative
    % 4. tap continues through the first point after that where the Z derivative goes positive
    % 5. tap stops at the first point after that where the Z derivative goes negative
    st = nanmedian(d.bvei(:,4));
    p0 = find(d.bvei(:,4) < st, 1);
    p1 = find(d.bvei(p0:end,4) > st+5, 1) + p0;
    p2 = find(diff(d.bvei(p1:end,4)) < 0, 1) + p1;
    p3 = find(diff(d.bvei(p2:end,4)) > 0, 1) + p2;
    p4 = find(diff(d.bvei(p3:end,4)) < 0, 1) + p3;
    hac = [p1 p2];
    imp = [p2 p4];
    mov = [1 p1];
    
    fprintf('mov=(%.1f,%.1f) imp=(%.1f,%.1f) ', mov(1)/3000, mov(2)/3000, imp(1)/3000, imp(2)/3000);
    
    wins = [.05 .25 1.25 6.25]*3000;
    d.romano = cell(1, length(wins));
    d.steinbach = cell(1, length(wins));

    for w=1:length(wins)
        fprintf('%d', wins(w));
        
        [pre, chunks] = romano_features('pre', d.biws, d.bvei, d.bai, mass, wins(w), 0, mov);
        d.romano{w} = pre;
        fprintf('r');

        % FIXME need wavelet toolbox for TR
        feats = {'MF', 'TR', 'SR', 'WV', 'SP', 'F', 'RG', 'Fr', 'FM'};
        d.steinbach{w} = zeros(size(chunks,1), length(feats) + 12*any(strcmp(feats, 'MF')));
        for i=1:size(chunks,1)
            mov_i = mov(1) + chunks(i,:);
            f = steinbach_features(feats, hac, imp, mov_i, struct('data', dft321(d.bai(:,2:4)), 'Fs', 3000), struct('data', sqrt(sum(d.biws(:,2:3).^2,2)), 'Fs', 3000), []);
            d.steinbach{w}(i,:) = f;
        end
        fprintf('s ');
    end
    
    data(materials{m}) = d;
    
    fprintf('done!\n');
end

%% analyze feature distributions

for m=1:length(materials)
    d = data(materials{m});
    d.rr2 = [];
    d.sr2 = [];
    d.mean_rr2 = [];
    d.mean_sr2 = [];
    
    for w=1:length(wins)
        r = d.romano{w};
        s = d.steinbach{w};
        rr2 = [];
        sr2 = [];
        for i=1:size(s,2)
            [~,~,~,~,stats] = regress(s(:,i), [r(:,[end-5 end-3]) ones(size(r,1),1)]);
            sr2(end+1) = stats(1);
        end
        for i=1:size(r,2)
            [~,~,~,~,stats] = regress(r(:,i), [r(:,[end-5 end-3]) ones(size(r,1),1)]);
            rr2(end+1) = stats(1);
        end
        d.rr2(end+1,:) = rr2;
        d.sr2(end+1,:) = sr2;
        d.mean_rr2(end+1) = nanmean(rr2);
        d.mean_sr2(end+1) = nanmean(sr2);
    end
    
    data(materials{m}) = d;
end

%% try ICRA stuff with steinbach features added

labels = {};
features = {};
split_idx = {};
for w=1:length(wins)
    labels{w} = [];
    features{w} = {{} []};
    for m=1:length(materials)
        d = data(materials{m});
        labels{w} = [labels{w}; ones(size(d.steinbach{w},1),1)*m];
        features{w}{1} = [features{w}{1}; d.romano{w}];
        features{w}{2} = [features{w}{2}; d.steinbach{w}];
    end
    
    split_idx{w} = randsample(1:2, length(labels{w}), true, [4/5 1/5]);
end

%%
    
nbins = 10;
binmode = 'naive';
alpha = 0.2;
nu = 0.1;
gamma = 1;

rconfusion = cell(length(wins), 1);
sconfusion = cell(length(wins), 1);
rsconfusion = cell(length(wins), 1);

for w=1:length(wins)
    train_features = [romano_features('post', features{w}{1}(split_idx{w}==1,:), nbins, binmode, alpha, 1) features{w}{2}(split_idx{w}==1,:)];
    train_labels = labels{w}(split_idx{w}==1,:);
    test_features = [romano_features('post', features{w}{1}(split_idx{w}==2,:), nbins, binmode, alpha, 1) features{w}{2}(split_idx{w}==2,:)];
    test_labels = labels{w}(split_idx{w}==2,:);
    
    trainmean = mean(train_features);
    trainrange = range(train_features);
    train_vectors = bsxfun(@rdivide, bsxfun(@minus, train_features, trainmean), trainrange);
    test_vectors = bsxfun(@rdivide, bsxfun(@minus, test_features, trainmean), trainrange);
    
    model = svmtrain(train_labels, train_vectors(:,1:end-21), sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', nu, gamma));
    answers = svmpredict(zeros(size(test_labels)), test_vectors(:,1:end-21), model, '-q');
    rconfusion{w} = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            rconfusion{w}(i,j) = nnz(answers(test_labels == i) == j);
        end
    end
    fprintf('Test set accuracy (win=%d, r): %g\n', wins(w), sum(diag(rconfusion{w}))/sum(sum(rconfusion{w})));

    model = svmtrain(train_labels, train_vectors(:,end-20:end), sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', nu, gamma));
    answers = svmpredict(zeros(size(test_labels)), test_vectors(:,end-20:end), model, '-q');
    sconfusion{w} = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            sconfusion{w}(i,j) = nnz(answers(test_labels == i) == j);
        end
    end
    fprintf('Test set accuracy (win=%d, s): %g\n', wins(w), sum(diag(sconfusion{w}))/sum(sum(sconfusion{w})));

    model = svmtrain(train_labels, train_vectors, sprintf('-m 1000 -s 1 -t 2 -n %g -g %g -q', nu, gamma));
    answers = svmpredict(zeros(size(test_labels)), test_vectors, model, '-q');
    rsconfusion{w} = zeros(length(materials));
    for i=1:length(materials)
        for j=1:length(materials)
            rsconfusion{w}(i,j) = nnz(answers(test_labels == i) == j);
        end
    end
    fprintf('Test set accuracy (win=%d, r+s): %g\n', wins(w), sum(diag(rsconfusion{w}))/sum(sum(rsconfusion{w})));
end

%% confusion matrix figures

rconf = figure;
imagesc(rconfusion{2});
colormap gray;
print -dpng rconf.png

sconf = figure;
imagesc(sconfusion{2});
colormap gray;
print -dpng sconf.png

%% print out results

fprintf('FEATURE SENSITIVITY TABLES\n');
fprintf('    \\begin{tabular}{c|ccccccccccc}\n');
fprintf('        Surface & Bin 1 & Bin 2 & Bin 3 & Bin 4 & Bin 5 & $\\overline{F_n}$ & $\\sigma(F_n)$ & $\\overline{v}$ & $\\sigma(v)$ & $\\overline{F_t}$ & $\\sigma(F_t)$ \\\\\n');
fprintf('        \\hline\n');
all_rr2 = [];
for m=1:length(materials)
    d = data(materials{m});
    
    fprintf('        %s', materials{m});
    for j=1:11
        fprintf(' & %.4f ', d.rr2(2,j));
    end
    fprintf('\\\\\n');
    all_rr2 = [all_rr2; d.rr2(2,:)];
end
fprintf('        \\emph{mean}');
for j=1:11
    fprintf(' & %.4f ', mean(all_rr2(:,j)));
end
fprintf('\\\\\n');
fprintf('    \\end{tabular}\n');

fprintf('    \\begin{tabular}{c|ccccccccc}\n');
fprintf('        Surface & AMCC & AMTR & AMSR & AMWV & AMSP & AMF & AMRG & Fr & FM \\\\\n');
fprintf('        \\hline\n');
all_sr2 = [];
for m=1:length(materials)
    d = data(materials{m});
    
    fprintf('        %s & %.4f', materials{m}, mean(d.sr2(2,1:13)));
    for j=14:21
        fprintf(' & %.4f ', d.sr2(2,j));
    end
    fprintf('\\\\\n');
    all_sr2 = [all_sr2; mean(d.sr2(2,1:13)) d.sr2(2,14:end)];
end
fprintf('        \\emph{mean}');
for j=1:9
    fprintf(' & %.4f ', mean(all_sr2(:,j)));
end
fprintf('\\\\\n');
fprintf('    \\end{tabular}\n');
