% CALC_CONV_MU_KERNEL Calculate convolution kernel for mu
%
% Usage
%    mu_kernel_f = calc_conv_mu_kernel_f(N, rot_matrices, ctfs, ctf_idx, ...
%        mu_est_opt);
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
%    mu_est_opt: A struct containing different options for the estimation of
%       the mean. The supported fields are:
%          - precision: The precision of the calculation. This is either
%            'double' or 'single' (default 'double').
%
% Output
%    mu_kernel_f: The Fourier transform of the convolution kernel
%       corresponding to the projection-backprojection operator of the
%       dataset. This is an array of size 2N-by-2N-by-2N.

function mu_kernel_f = calc_conv_mu_kernel_f(N, rot_matrices, ctfs, ...
       ctf_idx, mu_est_opt)
    if nargin < 5 || isempty(mu_est_opt)
        mu_est_opt = struct();
    end

    mu_est_opt = fill_struct(mu_est_opt, ...
        'precision', 'double');

mu_est_opt.precision

    sz_mult = 2;

    n = size(rot_matrices, 3);

    pts_rot = rotated_grids(N, rot_matrices);

    ctfsq = zeros([size(pts_rot, 2) size(pts_rot, 3)], mu_est_opt.precision);
    for s = 1:n
        ctfsq_s = ctfs{ctf_idx(s)}(frob_norm(pts_rot(:,:,s), 1)).^2;
		ctfsq_s = cast(ctfsq_s, mu_est_opt.precision);

		ctfsq(:,s) = ctfsq_s;
    end

    pts_rot = reshape(pts_rot, 3, []);
    ctfsq = reshape(ctfsq, 1, []);

    mu_kernel = (2/N)^2*anufft3(ctfsq, pi*pts_rot, ...
        sz_mult*N*ones(1,3));

    % Ensure symmetric kernel
    if mod(sz_mult*N, 2) == 0
        mu_kernel(1,:,:) = 0;
        mu_kernel(:,1,:) = 0;
        mu_kernel(:,:,1) = 0;
    end

    mu_kernel = fftn(ifftshift(ifftshift(ifftshift(mu_kernel, 1), 2), 3));
    mu_kernel_f = mu_kernel;
end
