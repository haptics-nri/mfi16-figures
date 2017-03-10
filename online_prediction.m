function [chunks, p, err, X, Y] = online_prediction(d, start, stop, varargin)

if length(varargin) == 1
    do_plot = false;
else
    do_plot = true;
    plot_title = varargin{1};
    legloc = varargin{2};
end

% unity-gain first-order low-pass filter
lpfa = .01;
lpfb = [1 lpfa-1];

force = [d.biws(start:stop,1) filtfilt(lpfa, lpfb, d.biws(start:stop,2:end))];
pos   = [d.bvei(start:stop,1) filtfilt(lpfa, lpfb, d.bvei(start:stop,2:end))];
speed = [d.bvei(start:stop,1) filtfilt(lpfa, lpfb, [diff(pos(:,2:end)) / mean(diff(pos(:,1))); zeros(1,size(pos,2)-1)])];

divs = start:301:stop; % make 100ms divisions

% extract chunks of data for each division
chunks = struct('t', [], 'in', [], 'out', []);
k = 0;
for i=1:length(divs)-1
    fprintf('constructing chunk %d/%d... \n', i+1, length(divs)-1);
    
    a = divs(i) - start + 1;
    b = divs(i+1)-1 - start + 1;
    
    t = force(a:b,1);
    Fn = force(a:b,4);
    Ft = sqrt(sum(force(a:b,2:3).^2,2));
    Vt = sqrt(sum(speed(a:b,2:3).^2,2));
    %Ft = zeros(size(Fn));
    %for j=1:length(Ft)
    %    Ft(j) = dot(force(a+j-1,2:3)', -speed(a+j-1,2:3)/sqrt(sum(speed(a+j-1,2:3).^2)));
    %end
    %Ft(Ft<0) = 0;
    
    if mean(Vt) >= 25 && max(Vt) <= 100
        k = k+1;
        chunks(k) = struct('t',   t, ...
                           'in',  [Fn Vt], ...
                           'out', Ft);
    end
end

% iteratively learn with increasing numbers of chunks, record all results
X = []; % training set input
Y = []; % training set output
p = []; % model parameters (one row per model)
err = []; % model error (one row per model, first column: error predicting the next chunk, second column: error predicting all chunks)
i = 0; while i < length(chunks)-1; i = i + 1;
    fprintf('predicting chunk %d/%d... \n', i+1, length(chunks));
    
    X = [X; chunks(i).in];
    Y = [Y; chunks(i).out];
    
    %b = robustfit(X, Y, 'logistic')'; % FIXME what kind of learner to use??
    b = (X\Y)';
    %[b, cs] = ransac_leastsquares(X, Y, 'sigma',0.25, 'verbose',false); b = b';
    p = [p; b];
    %err = [err; sqrt(mean((chunks(i+1).out - [ones(size(chunks(i+1).in,1), 1) chunks(i+1).in]*p(i,:)').^2))];
    err = [err; [sqrt(mean((chunks(i+1).out - chunks(i+1).in*p(i,:)').^2)) 0]];
    
    % throw out outliers
    if err(end,1) > 14
        fprintf('outlier\n');
        X = X(1:end-1,:);
        Y = Y(1:end-1,:);
        p = p(1:end-1,:);
        err = err(1:end-1,:);
        chunks(i+1) = [];
        i = i - 1;
    end
end
% go back and calculate overall error
for i=1:length(chunks)-1
    fprintf('evaluating chunk %d/%d... \n', i, length(chunks));
    %err(i,2) = sqrt(mean((reshape([chunks.out], [],1) - [ones(numel([chunks.out]), 1) cell2mat(arrayfun(@(s) s.in', chunks, 'uniformoutput', false))']*p(i,:)').^2));
    err(i,2) = sqrt(mean((reshape([chunks.out], [],1) - cell2mat(arrayfun(@(s) s.in', chunks, 'uniformoutput', false))'*p(i,:)').^2));
end

stopping = find(mean(movingmean(abs(diff(p)), 10), 2) < 1e-4, 1);

% plots
if do_plot
    subplot(4,1,1);
    plotyy(1:length(err), err(:,1), 1:length(err), err(:,2));
    legend('prediction', 'overall');
    xlabel('Iteration number');
    ylabel('Error', 'color','k');
    title(plot_title);
    subplot(4,1,2);
    plot(p);
    legend('normal force term', 'tangential speed term', 'location',legloc);
    xlabel('Iteration number');
    ylabel('Coefficient value');
    subplot(4,1,3);
    allin = cell2mat(arrayfun(@(s) s.in', chunks, 'uniformoutput', false))';
    plotyy(reshape([chunks.t], [],1), allin(:,1), ...
           reshape([chunks.t], [],1), allin(:,2));
    xlabel('Time (s)');
    ylabel('Learning inputs', 'color','k');
    legend('normal force (N)', 'tangential speed (mm/s)');
    subplot(4,1,4);
    plot(reshape([chunks.t], [],1), reshape([chunks.out], [],1), ...
         reshape([chunks.t], [],1), allin*p(10,:)', ...
         reshape([chunks.t], [],1), allin*p(stopping,:)', ...
         reshape([chunks.t], [],1), allin*p(end,:)');
    xlabel('Time (s)');
    ylabel('Tangential force (N)');
    legend('actual', 'predicted by 10th model', 'predicted at stopping point', 'predicted by final model')
end
