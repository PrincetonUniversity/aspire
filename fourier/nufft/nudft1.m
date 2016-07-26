% NUDFT1 Non-uniform discrete Fourier transform (1D)
%
% Usage
%    sig_f = nudft1(sig, fourier_pts);
%
% Input
%    sig: A vector of length containing a signal.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as a vector of size
%       K-by-1 with values in the range [-pi, pi].
%
% Output
%    sig_f: The Fourier transform of sig calculated at the specified freq-
%       uencies.
%
% See also
%    nufft1

function sig_f = nudft1(sig, fourier_pts)
	N = size(sig, 1);

	sig_f = zeros(size(fourier_pts, 1), 1);

	grid = ceil([-N/2:N/2-1]);

	for k = 1:size(fourier_pts, 1)
		sig_f(k) = exp(-i*(fourier_pts(k)*grid))*sig(:);
	end
end
