% FB_ROT_COVAR Calculate invariant covariance of Fourier-Bessel coefficients
%
% Usage
%    covar_coeff = fb_rot_covar(coeff, mean_coeff, basis, do_refl);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    mean_coeff: The coefficients of the mean image, as obtained from the
%       function fb_rot_mean.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%    do_refl: If true, enforce invariance to reflection (default false).
%
% Output
%    covar_coeff: The Fourier-Bessel coefficients of the covariance matrix in
%       the form of a basis.count-by-basis.count array. The covariance is
%       calculated from the images represented by the coeff array, along with
%       all possible rotations and reflections. As a result, the computed
%       covariance matrix is invariant to both reflection and rotation.

function covar_coeff = fb_rot_covar(coeff, mean_coeff, basis, do_refl)
    if nargin < 4 || isempty(do_refl)
        do_refl = true;
    end

    if basis.type ~= fb_basis_type() || ...
        numel(basis.sz) ~= 2

        error('Basis must be 2D Fourier-Bessel basis.');
    end

    covar_coeff = zeros(basis.count*ones(1, 2));

    ell = 0;
    mask = (basis.indices.ells == ell);

    coeff_ell = bsxfun(@minus, coeff(mask,:), mean_coeff(mask));

    covar_coeff(mask,mask) = coeff_ell*coeff_ell'/size(coeff, 2);

    for ell = 1:basis.ell_max
        mask = (basis.indices.ells == ell);
        mask_pos = mask & (basis.indices.sgns == +1);
        mask_neg = mask & (basis.indices.sgns == -1);

        covar_ell_diag = (coeff(mask_pos,:)*coeff(mask_pos,:)' + ...
            coeff(mask_neg,:)*coeff(mask_neg,:)')/(2*size(coeff, 2));

        covar_coeff(mask_pos,mask_pos) = covar_ell_diag;
        covar_coeff(mask_neg,mask_neg) = covar_ell_diag;

        if ~do_refl
            covar_ell_off = (coeff(mask_pos,:)*coeff(mask_neg,:)'/size(coeff, 2) - ...
               coeff(mask_neg,:)*coeff(mask_pos,:)')/(2*size(coeff, 2));

            covar_coeff(mask_pos,mask_neg) = covar_ell_off;
            covar_coeff(mask_neg,mask_pos) = covar_ell_off';
        end
    end
end
