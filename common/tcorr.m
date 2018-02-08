% TORR Tensor Frobenius correlation
%
% Usage
%    rho = tcorr(A, B, calc_dims);
%
% Input
%    A, B: Two arrays of the same size.
%    calc_dims: The dimensions along which to calculate the correlation
%       (default 1:ndims(A)).
%
% Output
%    rho: The correlation between A and B along the dimensions in calc_dims.

function rho = tcorr(A, B, calc_dims)
    if nargin < 3
        calc_dims = [];
    end

    sigma = tinner(A, B, calc_dims);
    rho = sigma./(tnorm(A, calc_dims).*tnorm(B, calc_dims));
end
