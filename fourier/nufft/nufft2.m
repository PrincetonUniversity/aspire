% NUFFT2 Wrapper for non-uniform FFT (2D)
%
% Usage
%    im_f = nufft2(im, fourier_pts);
%
% Input
%    im: An N-by-N array containing the pixel structure of an image.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as an N-by-2 array. These need to be in the
%       range [-pi, pi] in each dimension.
%
% Output
%    im_f: The Fourier transform of im at the frequencies fourier_pts.
%
% See also
%    nudft2

function im_f = nufft2(im, fourier_pts)
	p = nufft_initialize(size(im), size(fourier_pts, 1));

	p = nufft_set_points(p, fourier_pts);

	im_f = nufft_transform(p, im);

	nufft_finalize(p);
end
