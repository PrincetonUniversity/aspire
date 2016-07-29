% BROADCAST_IDX Generate broadcast indices
%
% Usage
%    idx = broadcast_idx(n, n_max);
%
% Input
%    n: The number of indices available.
%    n_max: The number of indices to broadcast over.
%
% Output
%    idx: A set of broadcast indices. If n is equal to n_max, this is simply
%       the range from 1 to n_max. Otherwise, it's a sequence of ones.

function idx = broadcast_idx(n, n_max)
    if n == n_max
        idx = 1:n_max;
    else
        idx = ones(1, n_max);
    end
end
