% REAL_SPH_HARMONIC Evaluate a real spherical harmonic
%
% Usage
%    y = real_sph_harmonic(j, m, theta, phi);
%
% Input
%    j, m: The order and degree of the spherical harmonic. These must satisfy
%        abs(m) < j.
%    theta, phi: The spherical coordinates of the points at which we want to
%       evaluate the real spherical harmonic. Here, `theta` is the latitude,
%       between 0 and pi, while `phi` is the longitude, between 0 and 2*pi.
%
% Output
%    y: The real spherical harmonics evaluated at the points (theta, phi).

% TODO: Specify conventions of these spherical harmonics in documentation.

function y = real_sph_harmonic(j, m, theta, phi)
	y = assoc_legendre(j, abs(m), cos(theta));

	y = sqrt((2*j+1)/(4*pi)/prod((j-abs(m)+1):(j+abs(m))))*y;

	if m < 0
		y = sqrt(2)*sin(abs(m)*phi).*y;
	elseif m > 0
		y = sqrt(2)*cos(abs(m)*phi).*y;
	end
end
