% BESSELJ_CLOSEST_ZERO Find closest zero using Halley's method
%
% Usage
%    z = besselj_closest_zero(nu, z0, max_iter);
%
% Input
%    nu: The order of the Bessel function.
%    z0: The initial value at which to start the search.
%    max_iter: Maximum number of iterations (default 10).
%
% Output
%    z: The zero of J_nu closest to z0 after max_iter iterations.

function z = besselj_closest_zero(nu, z0, max_iter)
	if nargin < 3 || isempty(max_iter)
		max_iter = 10;
	end

	z = z0;

	c = 8;

	for iter = 1:max_iter
		f = besselj(nu, z);
		if nu == 0
			fp = -besselj(1, z);
			fb = -1/2*(besselj(0, z)-besselj(2, z));
		elseif nu == 1
			fp = 1/2*(besselj(0, z)-besselj(2, z));
			fb = 1/2*(-besselj(1, z)-1/2*(besselj(1, z)-besselj(2, z)));
		else
			fp = 1/2*(besselj(nu-1, z)-besselj(nu+1, z));
			fb = 1/4*(besselj(nu-2, z)-besselj(nu+2, z));
		end

		dz = -2*f.*fp./(2*fp.^2-f.*fb);

		z = z + dz;

		if all(abs(dz) < c*eps(z))
			break;
		end
	end
end

