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
	if ndims(sig_f) > 2 || size(sig_f, 2) ~= 1
		error('Input ''sig_f'' must be of the form K-by-1.');
	end

	if ndims(fourier_pts) > 2 || any(size(fourier_pts) ~= [1 size(sig_f, 1)])
		error('Input ''fourier_pts'' must be of the form 1-by-K.');
	end

	if numel(sz) ~= 1 || any(floor(sz) ~= sz) || any(sz < 1)
		error('Input ''sz'' must be a positive integer scalar.');
	end

	N = sz(1);

	grid = ceil([-N/2:N/2-1]);

	sig = zeros(N, 1);

	for k = 1:size(grid, 2)
		sig(k) = exp(i*(grid(k)*fourier_pts))*sig_f(:);
	end
end
