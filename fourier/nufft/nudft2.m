% NUDFT2 Non-uniform discrete Fourier transform (2D)
%
% Usage
%    im_f = nudft2(im, fourier_pts);
%
% Input
%    im: An N-by-N array of pixels representing an image.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as a K-by-2 array, with
%       values in the range [-pi, pi].
%
% Output
%    im_f: The Fourier transform of im calculated at the specified freq-
%       uencies.
%
% See also
%    nufft2

function im_f = nudft2(im, fourier_pts)
	N = size(im, 1);

	if size(im, 2) ~= N
		error('only square images supported');
	end

	im_f = zeros(size(fourier_pts, 1), 1);

	grid = ceil([-N/2:N/2-1]);
	[grid_x, grid_y] = ndgrid(grid, grid);

	pts = [grid_x(:) grid_y(:)]';

	for k = 1:size(fourier_pts, 1)
		im_f(k) = exp(-i*(fourier_pts(k,:)*pts))*im(:);
	end
end
