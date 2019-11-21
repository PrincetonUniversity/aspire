function [ basis, sample_points, coeff_pos_k ] = Image_to_FB(projs)
% Compute the Fourier Bessel representation of the real space image data

% Input:
%   projs: subtomogram tilt slices for the i-th particle L*L*sN
% Output:
%   basis: 
%          basis.Phi_ns: bessel radial function J(R_{kq}\xi/c) for \xi in [0, c]
%          basis.ang_freqs: the associated angular frequencies (k in the paper)
%          basis.n_theta: the number of samples on concentric rings for polar Fourier transformation.
%   sample_points:
%          sample_points.r: positions in interval [0, c]
%          sample_points.w: weights
%   coeff_pos_k: Fourier bessel expansion coefficients in cell structure.
%   number of cells = k_max + 1. coeff_pos_k{i} contain coefficients for k
%   = i-1.
%
% Yuan Liu, 05/18/2016

% compute FB coefficients and functions
c = 0.5; % bandlimit (0, 0.5]
s = size(projs);
R = floor(size(projs,1)/2);  % Compact support radius in real space
num_pool = 4; % number of parallel pools
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_FB( n_r, R, c );
[ coeff_pos_k ]= FBcoeff_nfft(projs, R, basis, sample_points, num_pool);

end
