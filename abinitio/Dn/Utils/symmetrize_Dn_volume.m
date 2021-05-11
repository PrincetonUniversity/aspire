function dn_vol = symmetrize_Dn_volume(vol, symmetry_degree, order_2_gen, ...
    normalize_flag)
%SYMMETRIZE_DN_VOLUME Creates a Dn-symmetric volume from a given volume.
%   Detailed explanation goes here
%
%   Input:
%       vol - Some 3D volume.
%       symmetry_degree - an integer >= 3.
%       order_2_gen - Either the matrix g_x or g_y.
%       normalize_flag - A flag, weather the output volume is normalized or
%                        not.
%
%   Output:
%       dn_vol - A volume which approximately is a Dn-symmetric volume.
%       
%   Written by Elad Eatah April 2021.

cn_vol = symmetrize_Cn_volume(vol, symmetry_degree, normalize_flag);

if order_2_gen == diag([1, -1, -1])  % If gx
    xvol = fastrotate3x(cn_vol, 180);
elseif order_2_gen == diag([-1, 1, -1])  % If gy
    xvol = fastrotate3y(cn_vol, 180);
else
    error(strcat('The symmetry axis in the XY plane can ONLY',...
        ' be X-axis or Y-axis!'));
end

dn_vol = cn_vol + xvol;
if normalize_flag
    dn_vol = dn_vol / 2;
end
view3d(dn_vol);
end
