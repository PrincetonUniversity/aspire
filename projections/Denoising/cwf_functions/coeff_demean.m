% COEFF_DEMEAN Compute Fourier-Bessel coefficients for images
%
% Usage
%    coeff = coeff_demean(data, R, basis, sample_points, num_pool);
%
% Input
%    data: Array of size L-by-L-by-n containing images.
%    R: The radius of the disk on which the basis is supported.
%    basis: The basis struct, obtained from precomp_fb.
%    sample_points: The struct describing the quadrature, obtained from
%        precomp_fb.
%    num_pool: The number of parallel workers to use.
%
% Output
%    coeff: The Fourier-Bessel coefficients in the same format as given by
%       FBcoeff_nfft.

% Written by Tejal Bhamre based on Jane's FFB code - Oct 2015
% Reformatted, documented, and refactored by Joakim Anden - 2018-Apr-13

function coeff = coeff_demean(data, R, basis, sample_points, num_pool)
    n = size(data, 3);
    nb = floor(n/num_pool);
    remain = n - nb*num_pool;

    partition = [(nb+1)*ones(1, remain) nb*ones(1, num_pool-remain)];

    data = mat2cell(data, size(data, 1), size(data, 2), partition);

    coeff = FBcoeff_nfft(data, R, basis, sample_points);
end
