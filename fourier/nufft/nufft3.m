% NUFFT3 Wrapper for non-uniform FFT (3D)
%
% Usage
%    vol_f = nufft3(vol, fourier_pts);
%
% Input
%    vol: An N-by-N-by-N array containing the voxel structure of a volume.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 3-by-K array. These need to be in the
%       range [-pi, pi] in each dimension.
%
% Output
%    vol_f: The Fourier transform of vol at the frequencies fourier_pts.
%
% See also
%    nudft3

function vol_f = nufft3(vol, fourier_pts)
	p = nufft_initialize(size(vol), size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

	vol_f = nufft_transform(p, vol);

	nufft_finalize(p);
end
