function [Center,Radius,Error,CS] = sphereFit_ransac(X, sigma)

    if nargin == 1
        sigma = 0.075;
    end

    options = struct('sigma',      sigma, ...
                     'est_fun',    @est_fun, ...
                     'man_fun',    @man_fun, ...
                     'mode',       'RANSAC', ...
                     'reestimate', 'true');
    
    results = RANSAC(X', options);
    
    Center = results.Theta(1:3);
    Radius = results.Theta(4);
    Error  = results.Theta(5);
    CS     = results.CS;

end

function [theta, k] = est_fun(X, s)

    k = 4;

    if nargin == 0 || isempty(X)
        theta = [];
        return;
    end

    if nargin == 1 || isempty(s)
        s = 1:size(X,2);
    end
    
    [theta(1:3), theta(4), theta(5)] = sphereFit(X(:,s)');

end

% sum squared error for sphere
%
% sphere equation: (x - cx)^2 + (y - cy)^2 + (z - cz)^2 = r^2
%                  sum((xyz - c).^2, 2) = r^2
% E   = sum((xyz - c).^2, 2) - r^2
% SE  = (sum((xyz - c).^2, 2) - r^2).^2
% SSE = sum((sum((xyz - c).^2, 2) - r^2).^2)

function [E, T_noise_squared, d] = man_fun(theta, X, sigma, P_inlier)

    E = [];
    if ~isempty(theta) && ~isempty(X)
        E = (sum(bsxfun(@minus, X, theta(1:3)').^2) - theta(4)^2).^2;
            
        %fprintf('estimating error, E = %d x %d, sum(E) = %g\n', size(E,1), size(E,2), sum(E));
    end
    
    if nargout > 1
        if P_inlier == 0
            T_noise_squared = sigma;
        else
            d = 3;
            T_noise_squared = sigma^2 * chi2inv(P_inlier, d);
        end
    end

end
