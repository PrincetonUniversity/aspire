% FB_ROT_MEAN Calculate invariant mean of Fourier-Bessel coefficients
%
% Usage
%    mean_coeff = fb_rot_mean(coeff, basis);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%
% Output
%    mean_coeff: The mean coefficients of the images along with all possible
%       in-plane rotations and reflections. As a result, the mean image is
%       invariant to both reflection and rotation.

function mean_coeff = fb_rot_mean(coeff, basis)
    if basis.type ~= fb_basis_type() || ...
        numel(basis.sz) ~= 2

        error('Basis must be 2D Fourier-Bessel basis.');
    end

    mask = (basis.indices.ells == 0);

    mean_coeff = zeros(basis.count, 1);
    mean_coeff(mask,:) = mean(coeff(mask,:), 2);
end
