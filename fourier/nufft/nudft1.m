% NUDFT1 Non-uniform discrete Fourier transform (1D)
%
% Usage
%    sig_f = nudft1(sig, fourier_pts);
%
% Input
%    sig: An array of size N-by-1 containing a signal.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as a vector of size
%       1-by-K with values in the range [-pi, pi].
%
% Output
%    sig_f: The Fourier transform of sig calculated at the specified freq-
%       uencies.
%
% See also
%    nufft1

function sig_f = nudft1(sig, fourier_pts)
	if ndims(sig) > 2 || size(sig, 2) ~= 1
		error('Input ''sig'' must be of the form N-by-1.');
	end

	if ndims(fourier_pts) > 2 || size(fourier_pts, 1) ~= 1
		error('Input ''fourier_pts'' must be of the form 1-by-K.');
	end

	N = size(sig, 1);

	sig_f = zeros(size(fourier_pts, 2), 1);

	grid = ceil([-N/2:N/2-1]);

	for k = 1:size(fourier_pts, 2)
		sig_f(k) = exp(-i*(fourier_pts(k)*grid))*sig(:);
	end
end
