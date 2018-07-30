% TNORM Tensor Frobenius norm
%
% Usage
%    n = tnorm(A, calc_dims);
%
% Input
%    A: An array.
%    calc_dims: The dimensions along which to calculate the norm (default
%       1:ndims(A)).
%
% Output
%    n: The norm of A along the dimensions in calc_dims.

function n = tnorm(A, calc_dims)
    if nargin < 2
        calc_dims = [];
    end
    
    n = sqrt(tinner(A, A, calc_dims));
end
