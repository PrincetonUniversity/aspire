% ANUDFT3 Adjoint non-uniform discrete Fourier transform (3D)
%
% Usage
%    vol = anudft3(vol_f, fourier_pts, sz);
%
% Input
%    vol_f: A volume Fourier transform calculated at the frequencies specified
%       by fourier_pts. This is given as a vector.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are arranged as a 3-by-K array,
%       with values in the range [-pi, pi].
%    sz: The desired size of the output volume.
%
% Output
%    vol: The adjoint Fourier transform of vol_f at frequencies fourier_pts.

function vol = anudft3(vol_f, fourier_pts, sz)
	if ndims(vol_f) > 2 || size(vol_f, 2) ~= 1
		error('Input ''vol_f'' must be of the form K-by-1.');
	end

	if ndims(fourier_pts) > 2 || any(size(fourier_pts) ~= [3 size(vol_f, 1)])
		error('Input ''fourier_pts'' must be of the form 3-by-K.');
	end

	if numel(sz) ~= 3 || any(floor(sz) ~= sz) || any(sz < 1)
		error('Input ''sz'' must be a positive integer vector of length three.');
	end

	N = sz(1);

	if sz(2) ~= N || sz(3) ~= N
		error('Only cubic volumes supported.');
	end

	grid = ceil([-N/2:N/2-1]);
	[grid_x, grid_y, grid_z] = ndgrid(grid, grid, grid);

	pts = [grid_x(:) grid_y(:) grid_z(:)]';

	vol = zeros(N*ones(1, 3));

	for k = 1:size(pts, 2)
		vol(k) = exp(i*(pts(:,k)'*fourier_pts))*vol_f(:);
	end
end
