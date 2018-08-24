% FOURIER_BESSEL_BASIS Construct a Fourier-Bessel basis object (deprecated)
%
% Usage
%    basis = fourier_bessel_basis(sz, ell_max, domain);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only square images and cubic volumes are supported.
%    ell_max: The maximum order ell of the basis elements. If set to Inf,
%       the basis includes all ell such that the resulting basis vectors
%       are concentrated below the Nyquist frequency (default Inf).
%    domain: Specifies whether the decomposition should be in the spatial (0)
%       or frequency (1) domain (default 0).
%
% Output
%    basis: A Fourier-Bessel basis object corresponding to the parameters.
%
% Description
%    This function has been renamed to `fb_basis`. Please see the
%    documentation for that function for a description.

function basis = fourier_bessel_basis(varargin)
    warning('aspire:deprecated', ...
        ['`fourier_bessel_basis` is deprecated. Please call `fb_basis`.']);

    basis = fb_basis(varargin{:});
end
