% ANUFFT3 Wrapper for adjoint non-uniform FFT (3D)
%
% Usage
%    vol = anufft3(vol_f, fourier_pts, sz, nufft_opt);
%
% Input
%    vol_f: An array containing the Fourier transform of a volume at K
%       frequencies.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a K-by-3 array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting volume.
%    nufft_opt: A struct containing the fields:
%       - epsilon: The desired precision of the NUFFT (default 1e-15).
%
% Output
%    vol: The adjoint Fourier transform of vol_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft3

function vol = anufft3(vol_f, fourier_pts, sz, nufft_opt)
	if nargin < 4
		nufft_opt = [];
	end

	if adims(vol_f) < 1
		error('Input ''vol_f'' must be of the form K-by-L.');
	end

	if ndims(fourier_pts) > 2 || any(size(fourier_pts) ~= [3 size(vol_f, 1)])
		error('Input ''fourier_pts'' must be of the form 3-by-K.');
	end

	if numel(sz) ~= 3 || any(floor(sz) ~= sz) || any(sz < 1)
		error('Input ''sz'' must be a positive integer vector of length three.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 2), nufft_opt);

	p = nufft_set_points(p, fourier_pts);

	vol = nufft_adjoint(p, vol_f);

	nufft_finalize(p);
end
