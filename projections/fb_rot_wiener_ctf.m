% FB_ROT_WIENER_CTF Apply Wiener filter to filtered images
%
% Usage
%    coeff_est = fb_rot_wiener_ctf(coeff, filter_fb, filter_idx, ...
%       mean_coeff, covar_coeff, noise_var, basis);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    filter_fb: The block diagonal matrices representing filters arranged in a
%       cell array of length K. Each element is a block diagonal matrix
%       represented as another cell array that is manipulated by `blk_diag_*`
%       functions. Each block diagonal matrix corresponds to a basis.count-by-
%       basis-count size matrix.
%    filter_idx: The filter indices of the images in coeff as vector of length
%       n.
%    mean_coeff: The coefficients of the mean image, as obtained from the
%       function `fb_rot_mean_ctf`.
%    covar_coeff: A cell array representing the block diagonal covariance
%       matrix of the clean coefficients, as obtained from the function
%       `fb_rot_covar_ctf`.
%    noise_var: The variance of the noise in the images.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%
% Output
%    coeff_est: The estimated coefficients of the unfiltered images in the
%       Fourier-Bessel basis. These are obtained using a Wiener filter with
%       the specified covariance for the clean images and white noise of
%       variance `noise_var` for the noise.
%
% See also
%    fb_rot_mean_ctf, fb_rot_covar_ctf

function coeff_est = fb_rot_wiener_ctf(coeff, filter_fb, filter_idx, ...
    mean_coeff, covar_coeff, noise_var, basis)

    if ~ismember(basis.type, [fb_basis_type() ffb_basis_type()]) || ...
        numel(basis.sz) ~= 2

        error('Basis must be 2D Fourier-Bessel basis.');
    end

    blk_partition = blk_diag_partition(covar_coeff);
    precision = class(coeff);

    noise_covar_coeff = blk_diag_mult(noise_var, ...
        blk_diag_eye(blk_partition, precision));

    coeff_est = zeros(size(coeff), precision);

    for k = unique(filter_idx(:))'
        mask = (filter_idx == k);

        coeff_k = coeff(:,mask);

        filter_fb_k = filter_fb{k};
        filter_fb_k_t = blk_diag_transpose(filter_fb_k);

        sig_covar_coeff = ...
            blk_diag_mult(filter_fb_k, ...
            blk_diag_mult(covar_coeff, filter_fb_k_t));

        sig_noise_covar_coeff = blk_diag_add(sig_covar_coeff, ...
            noise_covar_coeff);

        mean_coeff_k = blk_diag_apply(filter_fb_k, mean_coeff);

        coeff_est_k = bsxfun(@minus, coeff_k, mean_coeff_k);
        coeff_est_k = blk_diag_solve(sig_noise_covar_coeff, coeff_est_k);
        coeff_est_k = blk_diag_apply( ...
            blk_diag_mult(covar_coeff, filter_fb_k_t), coeff_est_k);
        coeff_est_k = bsxfun(@plus, coeff_est_k, mean_coeff);

        coeff_est(:,mask) = coeff_est_k;
    end
end
