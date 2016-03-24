function [R,t,Error,CS] = rigid_ransac(A, B, sigma)

    if nargin == 2
        sigma = 0.1;
    end

    options = struct('sigma',      sigma, ...
                     'est_fun',    @est_fun, ...
                     'man_fun',    @man_fun, ...
                     'mode',       'RANSAC', ...
                     'reestimate', 'true');
    
    X = [A B];
    results = RANSAC(X', options);
    
    R      = reshape(results.Theta(1:9), [3 3]);
    t      = results.Theta(10:12)';
    Error  = sum(results.E);
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
    
    [R, t] = rigid_transform_3D(X(1:3,s)', X(4:6,s)');
    theta(1:9) = reshape(R, [1 9]);
    theta(10:12) = t;
end

function [E, T_noise_squared, d] = man_fun(theta, X, sigma, P_inlier)

    E = [];
    if ~isempty(theta) && ~isempty(X)
        R = reshape(theta(1:9), [3 3]);
        t = theta(10:12);
        E = sum((bsxfun(@plus, R\X(1:3,:), t') - X(4:6,:)).^2);
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
