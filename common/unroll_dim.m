% UNROLL_DIM Unroll trailing dimensions into one
%
% Usage
%    [X, roll_sz] = unroll_dim(X, dim);
%
% Input
%    X: The array whose trailing dimensions are to be unrolled.
%    dim: The dimension determining the start of the unrolling. All dimensions
%       after this will be unrolled into a single one.
%
% Output
%    X: The same signal X, but with the last ndims(X)-dim+1 dimensions rolled
%       into one.
%    roll_sz: The size of the unrolled dimensions.

function [X, roll_sz] = unroll_dim(X, dim)
    sz = [size(X) ones(1, max(0, dim-ndims(X)))];
    n_dim = length(sz);
    roll_sz = sz(dim:n_dim);

    X = reshape(X, [sz(1:dim-1) prod(roll_sz) 1]);
end
