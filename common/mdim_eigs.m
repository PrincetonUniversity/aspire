% MDIM_EIGS Multidimensional partial eigendecomposition
%
% Usage
%    [V, D] = mdim_eigs(A, k, sigma);
%
% Input
%    A: An array of size sig_sz-by-sig_sz, where sig_sz is a size containing
%       d dimensions. The array represents a matrix with d indices for its
%       rows and columns.
%    k: The number of eigenvalues and eigenvectors to calculate (default 6).
%    sigma: The type of eigenvalues to return ('lm', 'sm', 'la', 'sa', 'be',
%       'lr', 'sr', 'li', 'si', or a scalar). For more information, see the
%       documentation for eigs (default 'lm').
%
% Output
%    V: An array of eigenvectors of size sig_sz-by-k.
%    D: A matrix of size k-by-k containing the corresponding eigenvalues.

function [V, D] = mdim_eigs(A, k, sigma)
    if nargin < 2 || isempty(k)
        k = 6;
    end

    if nargin < 3 || isempty(sigma)
        sigma = 'lm';
    end

    if mod(ndims(A), 2) ~= 0
        error('Matrix must have an even number of dimensions');
    end

    d = ndims(A)/2;

    sig_sz = size(A);
    sig_sz = sig_sz(1:d);

    sig_len = prod(sig_sz);

    A = reshape(A, sig_len*ones(1, 2));

    [V, D] = eigs(A, k, sigma);

    V = reshape(V, [sig_sz k]);
end
