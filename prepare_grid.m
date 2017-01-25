% Prepares an index matrix for grid searching
% Each input argument is a 1D matrix or cell array
% Output is a prod(len(arg) for each arg)-x-nargin matrix where each row contains the indices to use for one grid search iteration
function gs_idx = prepare_grid(varargin)
    gs_limits = zeros(1,nargin);
    for i=1:nargin
        gs_limits(i) = length(varargin{i});
    end

    gs_idx = repmat(ones(size(gs_limits)), prod(gs_limits), 1);
    for i=2:size(gs_idx,1)
        gs_idx(i,:) = gs_idx(i-1,:);
        for j=size(gs_idx,2):-1:1
            if gs_idx(i,j) == gs_limits(j)
                gs_idx(i,j) = 1;
            else
                gs_idx(i,j) = gs_idx(i,j) + 1;
                break;
            end
        end
    end
end

