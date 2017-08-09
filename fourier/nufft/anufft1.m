% ANUFFT1 Wrapper for adjoint non-uniform FFT (1D)
%
% Usage
%    sig = anufft1(sig_f, fourier_pts, sz);
%
% Input
%    sig_f: An Fourier transform calculated at the K frequencies specified
%       by fourier_pts.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are in the form of an array size
%       1-by-K with values in the range [-pi, pi].
%    sz: The desired size of the output signal.
%
% Output
%    sig: The adjoint Fourier transform of sig_f at frequencies fourier_pts.

% See also
%    anudft1

function sig = anufft1(sig_f, fourier_pts, sz)
	if ndims(sig_f) > 2 || size(sig_f, 2) ~= 1
		error('Input ''sig_f'' must be of the form K-by-1.');
	end

	if ndims(fourier_pts) > 2 || any(size(fourier_pts) ~= [1 size(sig_f, 1)])
		error('Input ''fourier_pts'' must be of the form 1-by-K.');
	end

	if numel(sz) ~= 1 || any(floor(sz) ~= sz) || any(sz < 1)
		error('Input ''sz'' must be a positive integer scalar.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

	sig = nufft_adjoint(p, sig_f);

	nufft_finalize(p);
end
