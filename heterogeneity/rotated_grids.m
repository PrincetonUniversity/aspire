% ROTATED_GRIDS Generate rotated grids in 3D from rotation matrices
%
% Usage
%    [pts_rot, mask] = rotated_grids(N, rot_matrices, mask);
%
% Input
%    N: The resolution of the desired grids.
%    rot_matrices: An array of size 3-by-3-by-K containing K rotation matrices.
%    mask: A mask to apply to the grid. This has to be an N-by-N boolean array.
%       Alternatively, a scalar boolean can be used, in which case false gives
%       no mask, while true gives the standard disk mask provided by mesh_2d
%       (default true).
%    inclusive: Specifies whether the endpoints -1 and +1 should be included
%       in the rotated grids. For more information, see mesh_2d (default
%       false).
%
% Output
%    pts_rot: A set of rotated grids in three dimension as specified by the
%       rotation matrices and the mask.
%    mask: The mask used on the grids.

function [pts_rot, mask] = rotated_grids(N, rot_matrices, mask, inclusive)
    if nargin < 3
        mask = 1;
    end

    if nargin < 4 || isempty(inclusive)
        inclusive = false;
    end

    mesh2d = mesh_2d(N, inclusive);

    if numel(mask) == 1
        if mask == 0
            mask = true(N*ones(1,2));
        else
            mask = mesh2d.mask;
        end
    elseif size(mask) ~= N*ones(1,2)
        error('mask must be N by N');
    end

    num_pts = sum(mask(:));
    num_rots = size(rot_matrices, 3);

    pts = [mesh2d.x(mask) mesh2d.y(mask) zeros(num_pts, 1)]';
    pts_rot = zeros(3, num_pts, num_rots);
    for s = 1:num_rots
        pts_rot(:,:,s) = rot_matrices(:,:,s)*pts;
    end
end
