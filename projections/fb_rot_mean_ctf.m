% FB_ROT_MEAN_CTF Rotationally invariant mean from filtered images
%
% Usage
%    mean_coeff = fb_rot_mean_ctf(coeff, filter_fb, filter_idx, basis);
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
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%
% Output
%    mean_coeff: The mean coefficients of the images along with all possible
%       in-plane rotations and reflections. As a result, the mean image is
%       invariant to both reflection and rotation. The effect of the filters
%       in filter_fb are accounted for and inverted to yield a mean estimate
%       of the unfiltered images.
%
% See also
%    fb_rot_covar_ctf, fb_rot_mean, fb_rot_covar

function mean_coeff = fb_rot_mean_ctf(coeff, ctf_fb, ctf_idx, basis)
    if ~ismember(basis.type, [fb_basis_type() ffb_basis_type()]) || ...
        numel(basis.sz) ~= 2

        error('Basis must be 2D Fourier-Bessel basis.');
    end

    b = zeros(basis.count, 1);
    A = blk_diag_zeros(blk_diag_partition(ctf_fb{1}), ...
        class(ctf_fb{1}{1}));

    for k = unique(ctf_idx(:))'
        coeff_k = coeff(:,ctf_idx == k);

        weight = size(coeff_k, 2)/size(coeff, 2);

        mean_coeff_k = fb_rot_mean(coeff_k, basis);

        ctf_fb_k = ctf_fb{k};
        ctf_fb_k_t = blk_diag_transpose(ctf_fb_k);

        b = b + weight*blk_diag_apply(ctf_fb_k_t, mean_coeff_k(:));

        A = blk_diag_add(A, ...
            blk_diag_mult(weight, ...
            blk_diag_mult(ctf_fb_k_t, ctf_fb_k)));
    end

    mean_coeff = blk_diag_solve(A, b);
end
