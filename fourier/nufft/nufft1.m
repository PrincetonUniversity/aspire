% NUFFT1 Wrapper for non-uniform FFT (1D)
%
% Usage
%    sig_f = nufft1(sig, fourier_pts, nufft_opt);
%
% Input
%    sig: An array of size N-by-1 containing a signal.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as an array of size
%       1-by-K, with values in the range [-pi, pi].
%    nufft_opt: A struct containing the fields:
%       - epsilon: The desired precision of the NUFFT (default 1e-15).
%
% Output
%    sig_f: The Fourier transform of sig at the frequencies fourier_pts.
%
% See also
%    nudft1

function sig_f = nufft1(sig, fourier_pts, nufft_opt)
	if nargin < 3
		nufft_opt = [];
	end

	if adims(sig) < 1
		error('Input ''sig'' must be of the form N-by-L.');
	end

	if ndims(fourier_pts) > 2 || size(fourier_pts, 1) ~= 1
		error('Input ''fourier_pts'' must be of the form 1-by-K.');
	end

	p = nufft_initialize(asize(sig, 1), size(fourier_pts, 2), nufft_opt);

	p = nufft_set_points(p, fourier_pts);

	sig_f = nufft_transform(p, sig);

	nufft_finalize(p);
end
