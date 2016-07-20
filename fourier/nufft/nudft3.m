% NUDFT3 Non-uniform discrete Fourier transform (3D)
%
% Usage
%    vol_f = nudft3(vol, fourier_pts);
%
% Input
%    vol: An N-by-N-by-N array of voxels representing a volume.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as a 3-by-K array, with
%       values in the range [-pi, pi].
%
% Output
%    vol_f: The Fourier transform of vol calculated at the specified freq-
%       uencies.
%
% See also
%    nufft3

function vol_f = nudft3(vol, fourier_pts)
	N = size(vol, 1);

	if size(vol, 2) ~= N || size(vol, 3) ~= N
		error('only cube volumes supported');
	end

	vol_f = zeros(size(fourier_pts, 2), 1);

	grid = ceil([-N/2:N/2-1]);
	[grid_x, grid_y, grid_z] = ndgrid(grid, grid, grid);

	pts = [grid_x(:) grid_y(:) grid_z(:)]';

	for k = 1:size(fourier_pts, 2)
		vol_f(k) = exp(-i*(fourier_pts(:,k)'*pts))*vol(:);
	end
end
