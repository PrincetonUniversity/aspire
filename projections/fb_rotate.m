% FB_ROTATE Rotate image in Fourier-Bessel basis
%
% Usage
%    coeff_rot = fb_rotate(coeff, basis, angles);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%    angles: A vector of length n containing the angles of rotation between 0
%       and 2*pi.
%
% Output
%    coeff_rot: The coefficients of the original image, rotated by the
%       specified angles.

function coeff_rot = fb_rotate(coeff, basis, angles)
    coeff_rot = zeros(size(coeff, 1), max(size(coeff, 2), size(angles, 2)));

    mask = find(basis.indices.ells == 0);

    if size(coeff, 2) == 1 && size(angles, 2) > 1
        coeff_rot(mask,:) = repmat(coeff(mask), 1, size(angles, 2));
    else
        coeff_rot(mask,:) = coeff(mask,:);
    end

    for ell = 1:basis.ell_max
        is_ell = (basis.indices.ells == ell);

        mask_pos = find(is_ell & basis.indices.sgns == +1);
        mask_neg = find(is_ell & basis.indices.sgns == -1);

        cs = cos(ell*angles);
        si = sin(ell*angles);

        coeff_rot(mask_pos,:) = bsxfun(@times, cs, coeff(mask_pos,:)) - ...
            bsxfun(@times, si, coeff(mask_neg,:));
        coeff_rot(mask_neg,:) = bsxfun(@times, si, coeff(mask_pos,:)) + ...
            bsxfun(@times, cs, coeff(mask_neg,:));
    end
end
