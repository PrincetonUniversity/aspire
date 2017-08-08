% NUFFT3 Wrapper for non-uniform FFT (3D)
%
% Usage
%    vol_f = nufft3(vol, fourier_pts);
%
% Input
%    vol: An N1-by-N2-by-N3 array containing the voxel structure of a volume.
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
	if ndims(vol) ~= 3
		error('Input ''vol'' must be of the form N1-by-N2-by-N3.');
	end

	if ndims(fourier_pts) > 2 || size(fourier_pts, 1) ~= 3
		error('Input ''fourier_pts'' must be of the form 3-by-K.');
	end

	p = nufft_initialize(size(vol), size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

	vol_f = nufft_transform(p, vol);

	nufft_finalize(p);
end
