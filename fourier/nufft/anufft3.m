% ANUFFT3 Wrapper for adjoint non-uniform FFT (3D)
%
% Usage
%    vol = anufft3(vol_f, fourier_pts, sz);
%
% Input
%    vol_f: An array containing the Fourier transform of a volume at certain
%       points.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as an N-by-3 array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting volume.
%
% Output
%    vol: The adjoint Fourier transform of vol_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft3

function vol = anufft3(vol_f, fourier_pts, sz)
	if numel(sz) ~= 3
		error('Input ''sz'' must have three elements.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 1));

	p = nufft_set_points(p, fourier_pts);

	vol = nufft_adjoint(p, vol_f);

	nufft_finalize(p);
end
