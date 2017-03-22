% ANUDFT2 Adjoint non-uniform discrete Fourier transform (2D)
%
% Usage
%    im = anudft2(im_f, fourier_pts, sz);
%
% Input
%    im_f: An image Fourier transform calculated at the frequencies specified
%       by fourier_pts. This is given as a vector.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are arranged as a K-by-2 array,
%       with values in the range [-pi, pi].
%    sz: The desired size of the output image.
%
% Output
%    im: The adjoint Fourier transform of im_f at frequencies fourier_pts.

function im = anudft2(im_f, fourier_pts, sz)
	N = sz(1);

	if sz(2) ~= N
		error('only square images supported');
	end

	grid = ceil([-N/2:N/2-1]);
	[grid_x, grid_y] = ndgrid(grid, grid);

	pts = [grid_x(:) grid_y(:)]';

	im = zeros(N*ones(1, 2));

	for k = 1:size(pts, 2)
		im(k) = exp(i*(pts(:,k)'*fourier_pts'))*im_f(:);
	end
end
