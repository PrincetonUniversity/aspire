% SPH_BESSEL Compute spherical Bessel function values
%
% Usage
%    j = sph_bessel(ell, r)
%

function j = sph_bessel(ell, r)
	j = zeros(size(r));
	
	if ell == 0
		j(r==0) = 1;
	else
		j(r==0) = 0;
	end
	
	j(r~=0) = sqrt(pi./(2*r(r~=0))).*besselj(ell+1/2, r(r~=0));
end