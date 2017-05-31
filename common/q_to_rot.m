% Q_TO_ROT Convert quaternions to rotation matrices
%
% Usage
%    rot_matrices = q_to_rot(qs);
%
% Input
%    qs: An array of quaternions, arranged in columns. The size of this array
%       should be 4-by-n, where n is the number of quaternions.
%
% Output
%    rot_matrices: An array of size 3-by-3-by-n with the kth 3-by-3 block
%       equal to the rotation matrix corrsponding to the quaternion qs(:,k).

function rot_matrices = q_to_rot(qs)
	n = size(qs, 2);

	rot_matrices = zeros([3*ones(1, 2) n]);

	for k = 1:n
		q = qs(:,k);
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

		rot_matrices(:,:,k) = rot_matrix;
	end
end
