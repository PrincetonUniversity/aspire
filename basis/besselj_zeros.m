% BESSELJ_ZEROS Find k first zeros of J_nu
%
% Usage
%    z = besselj_zeros(nu, k);
%
% Input
%    nu: Order of the Bessel function J_nu. Must be less than 42.
%    k: The number of zeros desired.
%    max_iter: The maximum number of iterations used in Halley's method
%       for computing the exact zeros (default 10, see 
%       besselj_closest_zero).
%
% Output
%    z: The first k zeros of J_nu, ordered by increasing magnitude.
%
% Note
%    For now, the order nu has to be strictly less than 42 for the initial-
%    ization to work properly. Higher values give erroneous results and so are
%    not supported.

function z = besselj_zeros(nu, k, max_iter)
	if nu >= 42
		error('Only nu < 42 supported');
	end

	if nargin < 3
		max_iter = [];
	end

	% From J. McMahon, "On the Roots of the Bessel and Certain Related Func-
	% tions," Annals of Mathematics, Vol. 9, No. 1/6 (1894-1895), pp. 23-30.
	beta = 1/4*pi*(2*nu+4*[1:k])-1;
	m = 4*nu^2;
	z = beta;
	z = z-(m-1)./(8*beta)-4*(m-1)*(7*m-31)./(3*(8*beta).^3);
	z = z-32*(m-1)*(83*m^2-982*m+3779)./(15*(8*beta).^5);
	z = z-64*(m-1)*(6949*m^3-153855*m^2+1585743*m-6277237)./(105*(8*beta).^7);

	% Refine the initial estimates.
	z = besselj_closest_zero(nu, z, max_iter);
end

