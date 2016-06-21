%% weigh all sticks

datadir = '/Volumes/shared/Projects/Proton Pack/Data';
eps = {'20160620/weigh/1', '20160620/weigh/2', '20160617/weigh/1', '20160620/weigh/3'}; % TODO just scan the dir
sizes = [.25 .375 .5 .75];

masses = zeros(length(sizes),1);
coms = zeros(length(sizes),3);
bias = zeros(1,6);
for i=1:length(sizes)
    [masses(i), fbias, ~, coms(i,:), tbias, ~] = weigh({[datadir filesep eps{i}]});
    bias(1:3) = bias(1:3) + fbias;
    bias(4:6) = bias(4:6) + tbias';
end
bias = bias / length(sizes);

save stickweights masses coms bias

%% load texture data

date = '20160620';
epdirs = dir([datadir filesep date filesep 'stick']);
epdirs(arrayfun(@(e) e.name(1) == '.', epdirs)) = [];

eps = struct('endeff', cell(length(epdirs),1), ...
             'material', cell(length(epdirs),1), ...
             'flow', cell(length(epdirs),1), ...
             'data', cell(length(epdirs),1));
for i=1:length(epdirs)
    fprintf('[%d/%d] %s\n', i, length(epdirs), epdirs(i).name);
    j = str2double(epdirs(i).name);
    prefix = [datadir filesep date filesep 'stick' filesep epdirs(i).name];
    eps(j).flow = parse_flow([prefix filesep 'stick.flow']);
    eps(j).endeff = eps(j).flow.answers('tooling ball diameter').text;
    eps(j).material = eps(j).flow.answers('surface name').text;
    [v, f] = load_stick([prefix filesep]);
    eps(j).data = struct('vicon', v, 'force', f);
end

save stickdata eps

%% sync vicon/force data

for i=1:length(eps)
    fprintf('%d\n', i);
    %%
    v = eps(i).data.vicon;
    f = eps(i).data.force;
    
    % find the 4 taps in force (hand-tuned parameters seem to work)
    [fpks, flocs] = findpeaks(f(:,3), 'SortStr','descend', 'NPeaks',4, 'MinPeakProminence',2, 'MinPeakDistance',500);
    [flocs, idx] = sort(flocs); % resort in order
    fpks = fpks(idx);
    
    % find the 4 rises-before-the-taps in vicon
    [vpks, vlocs] = findpeaks(v(:,4), 'SortStr','descend', 'NPeaks',4, 'MinPeakDistance',10);
    [vlocs, idx] = sort(vlocs); % resort in order
    vpks = vpks(idx);
    
    eps(i).offset = mean(v(vlocs,1) - f(flocs,1));
end
