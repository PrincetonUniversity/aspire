% CALC_CONV_SIGMA_KERNEL Calculate convolution kernel for Sigma
%
% Usage
%    Sigma_kernel_f = calc_conv_Sigma_kernel_f(N, rot_matrices, ctfs, ...
%       ctf_idx, Sigma_est_opt);
%
% Input
%    N: The resolution of the images.
%    rot_matrices: A 3-by-3-by-n array of rotation matrices corresponding to
%       the viewing angles of the projections.
%    ctfs: Function handles of the CTFs in the dataset. This is a cell array,
%       with each element consisting of a function handle taking one input,
%       the radial frequencies at which the CTF is to be computed. For example
%       @(r)(r) represents the identity filter, where no CTF aberration is
%       present.
%    ctf_idx: The indices of the CTFs in the function handle array
%       corresponding to each projection. This is a vector of length n, where
%       n is the number of images.
%    Sigma_est_opt: A struct containing different options for the estimation of
%       the covariance. The supported fields are:
%          - precision: The precision of the calculation. This is either
%            'double' or 'single' (default 'double').
%
% Output
%    Sigma_kernel_f: The six-dimensional kernel corresponding to the
%       covariance projection-backprojection operator. This is an array of
%       size 2N-by-2N-by-2N-by-2N-by-2N-by-2N.

function Sigma_kernel_f = calc_conv_Sigma_kernel_f(N, rot_matrices, ctfs, ctf_idx, Sigma_est_opt)
    if nargin < 5
        Sigma_est_opt = struct();
    end

    factors = calc_conv_Sigma_kernel_f_factors(N, rot_matrices, ctfs, ctf_idx, Sigma_est_opt);

    factors = vol_to_vec(factors);

    Sigma_kernel_f = factors*factors';

    Sigma_kernel_f = vecmat_to_volmat(Sigma_kernel_f);
end
