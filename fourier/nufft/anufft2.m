% ANUFFT2 Wrapper for adjoint non-uniform FFT (2D)
%
% Usage
%    im = anufft2(im_f, fourier_pts, sz, nufft_opt);
%
% Input
%    im_f: An array containing the Fourier transform of an image at K
%       frequencies.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 2-by-K array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting image.
%    nufft_opt: A struct containing the fields:
%       - epsilon: The desired precision of the NUFFT (default 1e-15).
%
% Output
%    image: The adjoint Fourier transform of im_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft2

function im = anufft2(im_f, fourier_pts, sz, nufft_opt)
	if nargin < 4
		nufft_opt = [];
	end

	if adims(im_f) < 1
		error('Input ''im_f'' must be of the form K-by-L.');
	end

	if ndims(fourier_pts) > 2 || any(size(fourier_pts) ~= [2 size(im_f, 1)])
		error('Input ''fourier_pts'' must be of the form 2-by-K.');
	end

	if numel(sz) ~= 2 || any(floor(sz) ~= sz) || any(sz < 1)
		error('Input ''sz'' must be a positive integer vector of length two.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 2), nufft_opt);

	p = nufft_set_points(p, fourier_pts);

	im = nufft_adjoint(p, im_f);

	nufft_finalize(p);
end
