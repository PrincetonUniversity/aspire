% NUFFT1 Wrapper for non-uniform FFT (1D)
%
% Usage
%    sig_f = nufft1(sig, fourier_pts);
%
% Input
%    sig: A vector of length N containing a signal.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as an array of size 1-by-K,
%       with values in the range [-pi, pi].
%
% Output
%    sig_f: The Fourier transform of sig at the frequencies fourier_pts.
%
% See also
%    nudft1

function sig_f = nufft1(sig, fourier_pts)
	p = nufft_initialize(size(sig, 1), size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

	sig_f = nufft_transform(p, sig);

	nufft_finalize(p);
end
