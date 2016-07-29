% CALC_CONV_SIGMA_KERNEL_F_FACTORS Calculate convolution kernel factors for Sigma
%
% Usage
%    kernel_Sigma_f_factors = calc_conv_Sigma_kernel_f_factors( ...
%       N, rot_matrices, ctfs, ctf_idx, Sigma_est_opt);
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
%    kernel_Sigma_f_factors: The factors of the six-dimensional convolution
%       corresponding to the covariance projection-backprojection operator.
%       These are in an array of the form 2N-by-2N-by-2N-by-n, where n.

function factors = calc_conv_Sigma_kernel_f_factors(N, rot_matrices, ctfs, ctf_idx, Sigma_est_opt)
    global g_parallel;

    if nargin < 5 || isempty(Sigma_est_opt);
        Sigma_est_opt = struct();
    end

    Sigma_est_opt = fill_struct(Sigma_est_opt, ...
        'precision', 'double');

    n = size(rot_matrices, 3);

    sz_mult = 2;

    pts_rot = rotated_grids(N, rot_matrices);

    factors = zeros((sz_mult*N)^3, n, Sigma_est_opt.precision);

    if ~isempty(g_parallel) && g_parallel
        parfor s = 1:n
            factor_s = (2/N)^2*anufft3( ...
                ctfs{ctf_idx(s)}(frob_norm(pts_rot(:,:,s), 1)).^2, ...
                pi*pts_rot(:,:,s)', sz_mult*N*ones(1,3));

            % Ensure symmetric kernel
            if mod(sz_mult*N, 2) == 0
                factor_s(1,:,:) = 0;
                factor_s(:,1,:) = 0;
                factor_s(:,:,1) = 0;
            end

            factors(:,s) = 1/sqrt(n)*real(factor_s(:));
        end
    else
        for s = 1:n
            factor_s = (2/N)^2*anufft3( ...
                ctfs{ctf_idx(s)}(frob_norm(pts_rot(:,:,s), 1)).^2, ...
                pi*pts_rot(:,:,s)', sz_mult*N*ones(1,3));

            % Ensure symmetric kernel
            if mod(sz_mult*N, 2) == 0
                factor_s(1,:,:) = 0;
                factor_s(:,1,:) = 0;
                factor_s(:,:,1) = 0;
            end

            factors(:,s) = 1/sqrt(n)*real(factor_s(:));
        end
    end

    factors = vec_to_vol(factors);
    factors = ifftshift(ifftshift(ifftshift(factors, 1), 2), 3);
    factors = fft(fft(fft(factors, [], 1), [], 2), [], 3);
end
