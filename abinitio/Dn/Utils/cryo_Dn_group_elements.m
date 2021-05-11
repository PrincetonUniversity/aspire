function dn_elements = cryo_Dn_group_elements(symmetry_degree, order_2_gen)
%CRYO_DN_GROUP_ELEMENTS Generation of Dihedral group elements
%   This function generates the elements of the dihedral symmetry-group Dn.
%   Input:
%       symmetry_degree - an integer >= 2.
%       order_2_gen - The generator of order 2, a 3X3 matrix of rotation by
%                     Pi radians around some symmetry axis.
%
%   Output:
%       dn_elements - All group elements as 3X3X(2*symmetry_degree) matrix.
%       
%   Written by Elad Eatah March 2021. 

angles = linspace(0, 360, symmetry_degree + 1);
angles = angles(1:symmetry_degree);
cn_symmetries = quat2rotm(quaternion([zeros(numel(angles), 2), angles'], ...
    'eulerd', 'XYZ', 'frame'));
gx_symmetries = multiprod(order_2_gen, cn_symmetries, [1 2], [1 2]);
dn_elements = cat(3, cn_symmetries, gx_symmetries);
end
