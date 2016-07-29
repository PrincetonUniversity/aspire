% BACKPROJECT_CONV_SIGMA Backproject images for convolutional Sigma estimation
%
% Usage
%    [b_Sigma, b_Sigma_image, b_Sigma_noise] = backproject_conv_Sigma( ...
%        N, image_fun, rot_matrices, ctfs, ctf_idx, mu, sigma_noise, ...
%        Sigma_est_opt);
%
% Input
%    N: The resolution of the images.
%    image_fun: A function handle taking two inputs: s and t. Here, s is the
%       index of the first image and t is the number of images to extract.
%       The output is expect to be in a form of N-by-N-by-t.
%    rot_matrices: A 3-by-3-by-n array of rotation matrices corresponding to
%       the viewing angles of the projections given by image_fun.
%    ctfs: Function handles of the CTFs in the dataset. This is a cell array,
%       with each element consisting of a function handle taking one input,
%       the radial frequencies at which the CTF is to be computed. For example
%       @(r)(r) represents the identity filter, where no CTF aberration is
%       present.
%    ctf_idx: The indices of the CTFs in the function handle array
%       corresponding to each projection given by image_fun. This is a vector
%       of length n, where n is the number of images.
%    mu: The mean used to calculate the sample covariances of each image. This
%       can either be specified manually or estimated from the data.
%    sigma_noise: The standard deviation of the noise. Again, this can be
%       estimated from the data.
%    mu_est_opt: A struct containing different options for the estimation of
%       the covariance. The supported fields are:
%          - precision: The precision of the calculation. This is either
%            'double' or 'single' (default 'double').
%
% Output
%    b_Sigma: The backprojection of the sample covariances of the images,
%       with a noise term correction. This is an array of size N-by-N-by-N-
%       by-N-by-N-by-N.
%    b_Sigma_image: The image covariance part of b_Sigma.
%    b_Sigma_noise: The noise correction part of b_Sigma.

function [b_Sigma, b_Sigma_image, b_Sigma_noise] = ...
        backproject_conv_Sigma(N, image_fun, rot_matrices, ctfs, ctf_idx, ...
        mu, sigma_noise, Sigma_est_opt)

    if nargin < 8 || isempty(Sigma_est_opt)
        Sigma_est_opt = struct();
    end

    Sigma_est_opt = fill_struct(Sigma_est_opt, ...
        'precision', 'double');

    b_Sigma_image = backproject_conv_Sigma_image(N, image_fun, ...
        rot_matrices, ctfs, ctf_idx, mu, Sigma_est_opt);

    b_Sigma_noise = backproject_conv_Sigma_noise(N, rot_matrices, ctfs, ...
        ctf_idx, sigma_noise, Sigma_est_opt);

    b_Sigma = b_Sigma_image-b_Sigma_noise;
end

function b_Sigma_image = backproject_conv_Sigma_image(N, image_fun, ...
        rot_matrices, ctfs, ctf_idx, mu, Sigma_est_opt)
    global g_parallel;

    n = size(rot_matrices, 3);

    [pts_rot, mask] = rotated_grids(N, rot_matrices);

    factors = zeros(N^3, n, Sigma_est_opt.precision);

    im = image_fun(1, n);
    im_f = centered_fft2(im);

    if ~isempty(mu)
        proj_mu = vol_project(mu, rot_matrices);
        proj_mu = im_apply_radial_filter(proj_mu, ctfs(ctf_idx));

        proj_mu_f = centered_fft2(proj_mu);

        im_f = im_f-proj_mu_f;
    end

    if ~isempty(g_parallel) && g_parallel
        parfor s = 1:n
            im_s_f = im_f(:,:,s);
            im_s_f = im_s_f(mask);

            factor_s = 2/N*anufft3( ...
                im_s_f.*ctfs{ctf_idx(s)}( ...
                frob_norm(pts_rot(:,:,s), 1)), ...
                pi*pts_rot(:,:,s)', N*ones(1,3));

            factors(:,s) = real(factor_s(:));
        end
    else
        for s = 1:n
            im_s_f = im_f(:,:,s);
            im_s_f = im_s_f(mask);

            factor_s = 2/N*anufft3( ...
                im_s_f.*ctfs{ctf_idx(s)}( ...
                frob_norm(pts_rot(:,:,s), 1)), ...
                pi*pts_rot(:,:,s)', N*ones(1,3));

            factors(:,s) = real(factor_s(:));
        end
    end

    b_Sigma_image = factors*factors';

    b_Sigma_image = vecmat_to_volmat(b_Sigma_image);
end

function b_Sigma_noise = backproject_conv_Sigma_noise(N, rot_matrices, ...
        ctfs, ctf_idx, sigma_noise, Sigma_est_opt);
    n = size(rot_matrices, 3);

    pts_rot = rotated_grids(N, rot_matrices);

    ctfsq = zeros(size(pts_rot, 2), n, Sigma_est_opt.precision);
    for s = 1:n
        ctfsq(:,s) = ctfs{ctf_idx(s)}(frob_norm(pts_rot(:,:,s), 1)).^2;
    end

    ctfsq = bsxfun(@times, ctfsq, sigma_noise.^2);

    pts_rot = reshape(pts_rot, 3, []);
    ctfsq = reshape(ctfsq, 1, []);

    noise_kernel = (2/N)^2*anufft3(ctfsq, pi*pts_rot', ...
        2*N*ones(1,3));

    [I,J,K] = ind2sub(N*ones(1,3), 1:N^3);
    I = bsxfun(@minus, I', I)+N+1;
    J = bsxfun(@minus, J', J)+N+1;
    K = bsxfun(@minus, K', K)+N+1;
    ind = sub2ind((2*N)*ones(1,3), I(:), J(:), K(:));
    ind = reshape(ind, N^3*[1 1]);

    b_Sigma_noise = N^2*noise_kernel(ind);

    b_Sigma_noise = vecmat_to_volmat(b_Sigma_noise);
end
