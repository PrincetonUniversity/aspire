function rot_matrix = q_to_rot(q)
%
% Convert a quaternion into a rotation matrix.
%
%   Input: 
%           q: quaternion.
%   Output: 
%           rot_matrix: 3x3 rotation matrix

rot_matrix = zeros(3);

rot_matrix(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
rot_matrix(1,2) = 2*q(2)*q(3) - 2*q(1)*q(4); 
rot_matrix(1,3) = 2*q(1)*q(3) + 2*q(2)*q(4); 

rot_matrix(2,1) = 2*q(2)*q(3) + 2*q(1)*q(4);
rot_matrix(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2; 
rot_matrix(2,3) = -2*q(1)*q(2) + 2*q(3)*q(4); 

rot_matrix(3,1) = -2*q(1)*q(3) + 2*q(2)*q(4); 
rot_matrix(3,2) = 2*q(1)*q(2) + 2*q(3)*q(4); 
rot_matrix(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
end
