% FOURIER_BESSEL_BASIS_TYPE Type for the Fourier-Bessel basis (deprecated)
%
% Usage
%    type = fourier_bessel_basis_type();
%
% Output
%    type: The type identifier for the Fourier-Bessel basis.
%
% Description
%    This function has been renamed to `fb_basis_type`. Please see the
%    documentation for hat function for description.

function type = fourier_bessel_basis_type()
    warning('aspire:deprecated', ...
        ['`fourier_bessel_basis_type` is deprecated. Please call ' ...
         '`fb_basis_type`.']);

    type = fb_basis_type();
end
