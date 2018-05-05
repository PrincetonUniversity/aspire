% TINNER Tensor Frobenius inner product
%
% Usage
%    n = tinner(A, B, calc_dims);
%
% Input
%    A, B: Two arrays of the same size.
%    calc_dims: The dimensions along which to calculate the inner product
%       (default 1:ndims(A)).
%
% Output
%    n: The inner product between A and B along the dimensions in calc_dims.

function n = tinner(A, B, calc_dims)
    if nargin < 3
        calc_dims = [];
    end

    if ndims(A) ~= ndims(B) || any(size(A) ~= size(B))
        error('Input matrices must have the same shape.');
    end

    if isempty(calc_dims)
        calc_dims = 1:ndims(A);
    end

    sz = size(A);
    dims = 1:ndims(A);
    idx_dims = dims(~ismember(dims, calc_dims));

    A = permute(A, [calc_dims idx_dims]);
    B = permute(B, [calc_dims idx_dims]);
    
    A = reshape(A, prod(sz(calc_dims)), []);
    B = reshape(B, prod(sz(calc_dims)), []);

    n = sum(A.*conj(B), 1);

    n = reshape(n, [sz(idx_dims) 1 1]);
end
