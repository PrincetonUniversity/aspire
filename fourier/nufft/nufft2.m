% NUFFT2 Wrapper for non-uniform FFT (2D)
%
% Usage
%    im_f = nufft2(im, fourier_pts);
%
% Input
%    im: An N1-by-N2 array containing the pixel structure of an image.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 2-by-K array. These need to be in the
%       range [-pi, pi] in each dimension.
%
% Output
%    im_f: The Fourier transform of im at the frequencies fourier_pts.
%
% See also
%    nudft2

function im_f = nufft2(im, fourier_pts)
	if ndims(im) ~= 2
		error('Input ''im'' must be of the form N1-by-N2.');
	end

	if ndims(fourier_pts) > 2 || size(fourier_pts, 1) ~= 2
		error('Input ''fourier_pts'' must be of the form 2-by-K.');
	end

	p = nufft_initialize(size(im), size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

	im_f = nufft_transform(p, im);

	nufft_finalize(p);
end
