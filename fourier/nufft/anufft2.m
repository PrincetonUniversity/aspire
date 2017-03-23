% ANUFFT2 Wrapper for adjoint non-uniform FFT (2D)
%
% Usage
%    im = anufft2(im_f, fourier_pts, sz);
%
% Input
%    im_f: An array containing the Fourier transform of an image at certain
%       points.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as an N-by-2 array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting image.
%
% Output
%    image: The adjoint Fourier transform of im_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft2

function im = anufft2(im_f, fourier_pts, sz)
	if numel(sz) ~= 2
		error('Input ''sz'' must have two elements.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 1));

	p = nufft_set_points(p, fourier_pts);

	im = nufft_adjoint(p, im_f);

	nufft_finalize(p);
end
