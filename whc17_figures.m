%% generates figures for the WHC 2017 paper

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')
addpath('../steinbach/htk-mfcc');

DATADIR = '/home/haptics/shared/Projects/Proton Pack/Data';

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

date = '20170112';
load mini40_calib
[mass, fbias, ~, com, tbias] = weigh({fullfile(DATADIR, date, 'weigh', '1')});
tbias = tbias';
save mini40_calib fbias tbias;

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
flowtype = 'stickcam';

% load datasets from both collection days and merge them into one hashmap
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
    % 2. find first point where Z position goes 6mm above median (6mm threshold empirically determined)
    % 3. tap starts at the first point after that where the Z derivative goes negative
    % 4. tap continues through the first point after that where the Z derivative goes positive
    % 5. tap stops at the first point after that where the Z derivative goes negative
    st = nanmedian(d.bvei(:,4));
    p0 = find(d.bvei(:,4) < st, 1);
    p1 = find(d.bvei(p0:end,4) > st+6, 1) + p0;
    p2 = find(diff(d.bvei(p1:end,4)) < 0, 1) + p1;
    p3 = find(diff(d.bvei(p2:end,4)) > 0, 1) + p2;
    p4 = find(diff(d.bvei(p3:end,4)) < 0, 1) + p3;
    hac = [p1 p2];
    imp = [p2 p4];
    mov = [1 p1];
    
    fprintf('mov=(%.1f,%.1f) imp=(%.1f,%.1f) ', mov(1)/3000, mov(2)/3000, imp(1)/3000, imp(2)/3000);

    % calculate Romano + Steinbach features at three window lengths
    
    wins = [.05 .25 1.25]*3000;
    d.romano = cell(1, length(wins));
    d.steinbach = cell(1, length(wins));

    for w=1:length(wins)
        fprintf('%d', wins(w));
        
        % romano_features operates on the whole episode and choose chunks
        % it returns the chosen chunks so we can use exactly the same ones for steinbach_features
        [pre, chunks] = romano_features('pre', d.biws, d.bvei, d.bai, mass, wins(w), 0, mov);
        d.romano{w} = pre;
        fprintf('r');

        feats = {'MF', 'TR', 'SR', 'WV', 'SP', 'F', 'RG', 'Fr', 'FM'}; % all the AMxx features
        d.steinbach{w} = zeros(size(chunks,1), length(feats) + 12*any(strcmp(feats, 'MF')));
        % steinbach_features operates on one chunk at a time
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

%% train/test split

% Create 80% train/20% test splits for each window length.
% outputs: labels{}, features{}, cv{}

labels = {};
features = {};
cv = {};
for w=1:length(wins)
    labels{w} = [];
    features{w} = {{} []};
    for m=1:length(materials)
        d = data(materials{m});
        labels{w} = [labels{w}; ones(size(d.steinbach{w},1),1)*m];
        features{w}{1} = [features{w}{1}; d.romano{w}];
        features{w}{2} = [features{w}{2}; d.steinbach{w}];
    end
    
    cv{w} = cvpartition(length(labels{w}), 'KFold', 5);
end

%% grid searches

% grid search at every window length to find the best hyperparameters for Romano-only, Steinbach-only, and combined
% gs_acc{w} is a KxNx3 matrix: K = outer-CV splits, N = hyperparam combos, 3 = R, S, R+S, values = mean inner-CV accuracy
    
nbins = [3 5 10];
binmode = {'naive' 'perceptual'};
alpha = [0.2 0.3 0.4];
nu = [0.05 0.1 0.15];
gamma = [1 10];

gs_acc = cell(1,length(wins));
for w=1:length(wins)
    gs_acc{w} = whc17_grid(materials, labels{w}, features{w}, cv{w}, nbins, binmode, alpha, nu, gamma);
end

%% analyze feature distributions

% Analyzes the correlation between scan-time parameters and each Romano/Steinbach feature.
% The scan-time parameters are mean normal force and mean tip speed, which are at end-5 and end-3 of the Romano feature vector.
% For analysis, we recalculate the Romano features with the best parameters found via grid search.

for m=1:length(materials)
    d = data(materials{m});
    d.rr2 = [];
    d.sr2 = [];
    d.mean_rr2 = [];
    d.mean_sr2 = [];
    
    for w=1:length(wins)
        % recalculate using grid search results

        % use the mode over CV folds of the best grid indices
        gsi = zeros(1,cv{w}.NumTestSets);
        for i=1:cv{w}.NumTestSets
            [~,gsi] = max(gs_acc{w}(i,:,1));
        end
        gsi = mode(gsi);
        gp = grid_params(gsi,:);
        r = romano_features('post', d.romano{w}, nbins(gp(1)), binmode{gp(2)}, alpha(gp(3)), 1);
        s = d.steinbach{w};

        % calculate correlation of each feature
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

        % store back to hashmap
        d.rr2(end+1,:) = rr2;
        d.sr2(end+1,:) = sr2;
        d.mean_rr2(end+1) = nanmean(rr2);
        d.mean_sr2(end+1) = nanmean(sr2);
    end
    
    data(materials{m}) = d;
end

%% evaluate best models and make figures

rconfusion = cell(1,length(wins));
sconfusion = cell(1,length(wins));
rsconfusion = cell(1,length(wins));
racc = cell(1,length(wins));
sacc = cell(1,length(wins));
rsacc = cell(1,length(wins));

grid_params = prepare_grid(nbins, binmode, alpha, nu, gamma);

for w=1:length(wins)
    [rconfusion{w}, racc{w}] = whc17_test(materials, wins(w), features{w}, labels{w}, cv{w}, 'romano', {nbins binmode alpha nu gamma}, grid_params, gs_acc{w}(:,:,1), sprintf('rconf%d.pdf', w));

    [sconfusion{w}, sacc{w}] = whc17_test(materials, wins(w), features{w}, labels{w}, cv{w}, 'steinbach', {nbins binmode alpha nu gamma}, grid_params, gs_acc{w}(:,:,2), sprintf('sconf%d.pdf', w));

    [rsconfusion{w}, rsacc{w}] = whc17_test(materials, wins(w), features{w}, labels{w}, cv{w}, 'both', {nbins binmode alpha nu gamma}, grid_params, gs_acc{w}(:,:,3), sprintf('rsconf%d.pdf', w));
end

%% construct latex table from grid search results

fprintf('\n\nTABLE 1\n\n');
indent = '            ';
fprintf('%s\\begin{tabular}{c|cccc}\n', indent);
fprintf('%s    Window & Scan-dependent (\\%%) & Scan-free (\\%%) & All \\\\\n', indent);
fprintf('%s    \\hline\n', indent);
for w=1:length(wins)
    elem = '$%2.2f~(\\sigma = %1.1f)$';
    fprintf(['%s    \\units{%1.2f}{s} & ' elem ' & ' elem ' & ' elem], ...
        indent, wins(w)/3000, ...
        mean(racc{w})*100, std(racc{w})*100, ...
        mean(sacc{w})*100, std(sacc{w})*100, ...
        mean(rsacc{w})*100, std(rsacc{w})*100);
    if w ~= length(wins)
        fprintf(' \\\\\n');
    else
        fprintf('\n');
    end
end
fprintf('%s\\end{tabular}\n', indent);

%% feature correlation figures

% visualizes the results of the feature correlation calculations using imagesc

% collect values out of the hashmap
featsens = cell(length(wins),2);
for w=1:3
    featsens{w,1} = zeros(length(materials), 16);
    featsens{w,2} = zeros(length(materials), 9);

    for m=1:length(materials)
        d = data(materials{m});
        featsens{w,1}(m,:) = d.rr2(w,:);
    end
    for m=1:length(materials)
        d = data(materials{m});
        featsens{w,2}(m,:) = [mean(d.sr2(w,1:13)) d.sr2(w,14:end)];
    end
end

% make figures
for w=1:length(wins)
    % ROMANO FIGURE
    imagesc([featsens{w,1}; mean(featsens{w,1})]);
    colormap(flipud(gray));
    ax = gca;
    ax.XTick = 1:size(featsens{w,1},2);
    ax.XTickLabel = [];
    xticklabels = [arrayfun(@(n) sprintf('Bin %d', n), 1:(size(featsens{w,1},2)-6), 'uniformoutput',false) {'$\overline{F_n}$' '$\sigma(F_n)$' '$\overline{v}$' '$\sigma(v)$' '$\overline{F_t}$' '$\sigma(F_t)$'}];
    for i=1:length(xticklabels)
        text(i, ax.YLim(2) + 0.5, xticklabels{i}, 'interpreter','latex', 'horizontalalignment','right', 'rotation',90);
    end
    ax.XTickLabelRotation = 90;
    ax.YTick = 1:(size(featsens{w,1},1)+1);
    ax.YTickLabel = [materials 'mean'];
    ax.FontName = 'Courier';
    xl = xlabel('Feature', 'fontname','arial', 'fontsize',16);
    xly = xl.Position(2) + 3;
    xl.Position(2) = xly;
    ylabel('Surface', 'fontname','arial', 'fontsize',16);
    print('-dpdf', sprintf('rsens%d.pdf', w));

    % STEINBACH FIGURE
    imagesc([featsens{w,2}; mean(featsens{w,2})]);
    colormap(flipud(gray));
    ax = gca;
    ax.XTick = 1:size(featsens{w,2},2);
    ax.XTickLabel = {'AMCC' 'AMTR' 'AMSR' 'AMWV' 'AMSP' 'AMF' 'AMRG' 'Fr' 'FM'};
    ax.XTickLabelRotation = 90;
    ax.YTick = 1:(size(featsens{w,2},1)+1);
    ax.YTickLabel = [materials 'mean'];
    ax.FontName = 'Courier';
    xl = xlabel('Feature', 'fontname','arial', 'fontsize',16);
    xl.Position(2) = xly;
    ylabel('Surface', 'fontname','arial', 'fontsize',16);
    print('-dpdf', sprintf('ssens%d.pdf', w));
end

%% replicate steinbach results

steinbach_classifier

