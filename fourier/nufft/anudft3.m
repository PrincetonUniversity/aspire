% ANUDFT3 Adjoint non-uniform discrete Fourier transform (3D)
%
% Usage
%    vol = anudft3(vol_f, fourier_pts);
%
% Input
%    vol_f: A volume Fourier transform calculated at the frequencies specified
%       by fourier_pts. This is given as a vector.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are arranged as a K-by-3 array,
%       with values in the range [-pi, pi].
%
% Output
%    vol: The adjoint Fourier transform of vol_f at frequencies fourier_pts.

function vol = anudft3(vol_f, fourier_pts, sz)
	N = sz(1);

	if sz(2) ~= N || sz(3) ~= N
		error('only cube volumes supported');
	end

	grid = ceil([-N/2:N/2-1]);
	[grid_x, grid_y, grid_z] = ndgrid(grid, grid, grid);

	pts = [grid_x(:) grid_y(:) grid_z(:)]';

	vol = zeros(N*ones(1, 3));

	for k = 1:size(pts, 2)
		vol(k) = exp(i*(pts(:,k)'*fourier_pts'))*vol_f(:);
	end
end
