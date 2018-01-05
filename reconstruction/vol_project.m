% VOL_PROJECT Project volume along rotations
%
% Usage
%    im = vol_project(vol, rot_matrices, half_pixel);
%
% Input
%    vol: An L-by-L-by-L array containing the voxel structure of a volume.
%    rot_matrices: A set of rotation matrices of the form 3-by-3-by-n, corre-
%       sponding to n different projection directions.
%    half_pixel: If true, centers the rotation around a half-pixel (default false).
%
% Output
%    im: An L-by-L-by-n array containing the projections of the volumes in the
%       specified directions.

function im = vol_project(vol, rot_matrices, half_pixel)
    if nargin < 3 || isempty(half_pixel)
        half_pixel = false;
    end

    L = size(vol, 1);

    if numel(size(vol)) ~= 3
        error('Input `vol` must have three dimensions.');
    end

    if size(vol, 2) ~= L || size(vol, 3) ~= L
        error(['Input `vol` must be an array of the form ' ...
            'L-by-L-by-L.']);
    end

    if mod(L, 2) == 0 && half_pixel
        [X, Y] = ndgrid(-L/2:L/2-1, -L/2:L/2-1);
    end

    n = size(rot_matrices, 3);

    pts_rot = rotated_grids(L, rot_matrices, half_pixel);

    pts_rot = reshape(pts_rot, [3 L^2*n]);

    im_f = nufft3(vol, pts_rot);

    im_f = reshape(im_f, [L*ones(1, 2) n]);

    if mod(L, 2) == 0
        if ~half_pixel
            im_f(1,:,:) = 0;
            im_f(:,1,:) = 0;
        else
            phase_shift = -reshape(sum(pts_rot, 1), [L*ones(1, 2) n])/2;
            phase_shift = phase_shift + 2*pi*(X+Y+1)/(2*L);
            im_f = bsxfun(@times, im_f, exp(1i*phase_shift));
        end
    end

    im = icfft2(im_f);

    if mod(L, 2) == 0 && half_pixel
        phase_shift = 2*pi*(X+Y)/(2*L);
        im = bsxfun(@times, im, exp(1i*phase_shift));
    end

    im = real(im);

    im = permute(im, [2 1 3]);
end
