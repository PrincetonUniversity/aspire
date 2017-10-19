% ROTATED_GRIDS Generate rotated Fourier grids in 3D from rotation matrices
%
% Usage
%    pts_rot = rotated_grids(N, rot_matrices, half_pixel);
%
% Input
%    N: The resolution of the desired grids.
%    rot_matrices: An array of size 3-by-3-by-K containing K rotation matrices.
%    half_pixel: If true, centers the rotation around a half-pixel (default false).
%
% Output
%    pts_rot: A set of rotated Fourier grids in three dimensions as specified by
%       the rotation matrices. Frequencies are in the range [-pi, pi].

function pts_rot = rotated_grids(L, rot_matrices, half_pixel)
    if nargin < 3 || isempty(half_pixel)
        half_pixel = false;
    end

    mesh2d = mesh_2d(L);

    if half_pixel
        mesh2d.x = mesh2d.x + 1/L;
        mesh2d.y = mesh2d.y + 1/L;
    end

    num_pts = L^2;
    num_rots = size(rot_matrices, 3);

    pts = pi*[mesh2d.x(:) mesh2d.y(:) zeros(num_pts, 1)]';
    pts_rot = zeros(3, num_pts, num_rots);
    for s = 1:num_rots
        pts_rot(:,:,s) = rot_matrices(:,:,s)*pts;
    end

    pts_rot = reshape(pts_rot, [3 L*ones(1, 2) num_rots]);
end
