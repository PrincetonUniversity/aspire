% FB_REFLECT Reflect image along x axis in Fourier-Bessel basis
%
% Usage
%    coeff_refl = fb_reflect(coeff, basis, refl);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%    refl: A vector of length n specifying which images should be reflected.
%       A value of zero means no reflection while a one means reflection.
%
% Output
%    coeff_refl: The coefficients of the original images, reflected along the
%       x axis (the first dimension).

function coeff_refl = fb_reflect(coeff, basis, refl)
    mask_neg = (basis.indices.sgns == -1);

    mult = (-1).^((basis.indices.ells + double(mask_neg))*refl);

    coeff_refl = bsxfun(@times, coeff, mult);
end
