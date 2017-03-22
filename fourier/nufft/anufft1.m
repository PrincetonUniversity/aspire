% ANUFFT1 Wrapper for adjoint non-uniform FFT (1D)
%
% Usage
%    sig = anufft1(sig_f, fourier_pts, sz);
%
% Input
%    sig_f: An Fourier transform calculated at the frequencies specified
%       by fourier_pts.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are in the form of a vector of
%       length K with values in the range [-pi, pi].
%    sz: The desired size of the output signal.
%
% Output
%    sig: The adjoint Fourier transform of sig_f at frequencies fourier_pts.

% See also
%    anudft1

function sig = anufft1(sig_f, fourier_pts, sz)
	p = nufft_initialize(sz, size(fourier_pts, 1));

	p = nufft_set_points(p, fourier_pts);

	sig = nufft_adjoint(p, sig_f);

	nufft_finalize(p);
end
