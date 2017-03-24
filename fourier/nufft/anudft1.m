% ANUDFT1 Adjoint non-uniform discrete Fourier transform (1D)
%
% Usage
%    sig = anudft1(sig_f, fourier_pts, sz);
%
% Input
%    sig_f: A Fourier transform calculated at the K frequencies specified
%       by fourier_pts. Must be an array of size K-by-1.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are in the form of a vector of
%       size 1-by-K with values in the range [-pi, pi].
%    sz: The desired size of the output signal.
%
% Output
%    sig: The adjoint Fourier transform of sig_f at frequencies fourier_pts.

function sig = anudft1(sig_f, fourier_pts, sz)
	N = sz(1);

	grid = ceil([-N/2:N/2-1]);

	sig = zeros(N, 1);

	for k = 1:size(grid, 2)
		sig(k) = exp(i*(grid(k)*fourier_pts))*sig_f(:);
	end
end
