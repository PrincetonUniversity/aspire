% BACKPROJECT_CONV_MU Backproject images for convolutional mu estimation
%
% Usage
%    b_mu = backproject_conv_mu(N, image_fun, rot_matrices, ctfs, ctf_idx, ...
%       mu_est_opt);
%
% Input
%    N: The resolution of the images.
%    image_fun: A function handle taking two inputs: s and t. Where s is the
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
%    mu_est_opt: A struct containing different options for the estimation of
%       the mean. The supported fields are:
%          - precision: The precision of the calculation. This is either
%            'double' or 'single' (default 'double').
%
% Output
%    b_mu: The backprojection of the images in volume space. This is a cubic
%       volume array of size N-by-N-by-N.

function b_mu = backproject_conv_mu(N, image_fun, rot_matrices, ctfs, ...
       ctf_idx, mu_est_opt)
    if nargin < 6 || isempty(mu_est_opt)
        mu_est_opt = struct();
    end

    mu_est_opt = fill_struct(mu_est_opt, ...
        'precision', 'double');

    n = size(rot_matrices, 3);

    [pts_rot, mask] = rotated_grids(N, rot_matrices);

    im = image_fun(1, n);

    im = cast(im, mu_est_opt.precision);

    im = im_apply_radial_filter(im, ctfs(ctf_idx));
    im_f = im_to_vec(centered_fft2(im));
    im_f = im_f(mask,:);

    pts_rot = reshape(pts_rot, 3, []);
    im_f = reshape(im_f, 1, []);

    b_mu = 1/n*2/N*anufft3(im_f, pi*pts_rot', N*ones(1,3));
end
