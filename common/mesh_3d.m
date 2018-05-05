% MESH_3D Define a volume mesh and mask
%
% Usage
%    mesh = mesh_3d(L, inclusive);
%
% Input
%    L: The dimension of the desired mesh.
%    inclusive: Specifies whether both endpoints -1 and +1 should be included
%       in each axis. If true, -1 and +1 are always included, and so 0 is only
%       included for odd L. If true, -1 is included for even L, while for odd
%       L, neither -1 or +1 is included (default false).
%
% Output
%    mesh: A structure containing the fields:
%          - x, y and z: x, y and z coordinates in an L-by-L-by-L grid.
%          - r, theta, phi: Spherical coordinates of the points on the grid.

function mesh = mesh_3d(L, inclusive)
    if nargin < 2 || isempty(inclusive)
        inclusive = false;
    end

    persistent p_meshes;

    if size(p_meshes, 1) >= L && ...
        size(p_meshes, 2) >= int32(inclusive)+1 && ...
        ~isempty(p_meshes{L, int32(inclusive)+1})
        mesh = p_meshes{L, int32(inclusive)+1};
        return;
    end

    if ~inclusive
        grid = ceil([-L/2:L/2-1])/(L/2);
    else
        grid = [-(L-1)/2:(L-1)/2]/((L-1)/2);
    end

    [mesh.x,mesh.y,mesh.z] = ndgrid(grid, grid, grid);

    [mesh.phi, mesh.theta, mesh.r] = cart2sph(mesh.x, mesh.y, mesh.z);

    mesh.theta = pi/2-mesh.theta;

    p_meshes{L, int32(inclusive)+1} = mesh;
end
