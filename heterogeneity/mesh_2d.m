% MESH_2D Define an image mesh and mask
%
% Usage
%    mesh = mesh_2d(N, inclusive);
%
% Input
%    N: The dimension of the desired mesh.
%    inclusive: Specifies whether both endpoints -1 and +1 should be included
%       in each axis. If true, -1 and +1 are always included, and so 0 is only
%       included for odd N. If true, -1 is included for even N, while for odd
%       N, neither -1 or +1 is included (default false).
%
% Output
%    mesh: A structure containing the fields:
%          - x and y: x and y coordinates in an N-by-N grid.
%          - mask: A disk mask containing only the central disk of the grid.
%          - r, phi: Polar coordinates of the points in the disk mask.

function mesh = mesh_2d(N, inclusive)
    if nargin < 2 || isempty(inclusive)
        inclusive = false;
    end

    persistent p_meshes;

    if size(p_meshes, 1) >= N && ...
        size(p_meshes, 2) >= int32(inclusive)+1 && ...
        ~isempty(p_meshes{N, int32(inclusive)+1})
        mesh = p_meshes{N, int32(inclusive)+1};
        return;
    end

    if ~inclusive
        grid = ceil([-N/2:N/2-1])/(N/2);
    else
        grid = [-(N-1)/2:(N-1)/2]/((N-1)/2);
    end
    % NOTE: Need to use ndgrid because meshgrid swaps x and y...
    [mesh.x,mesh.y] = ndgrid(grid, grid);

    mesh.mask = (mesh.x.^2+mesh.y.^2)<1-1e-10;

    [mesh.phi, mesh.r] = cart2pol(mesh.x(mesh.mask), mesh.y(mesh.mask));

    p_meshes{N, int32(inclusive)+1} = mesh;
end
