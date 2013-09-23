function [inv_rot_matrices]=angle_to_invrot(angles)
% This function convert the Euler angles to inverse rotation matrices.
n_proj=size(angles,3);
rot_matrices = zeros(3,3,n_proj);

t=zeros(3);t(1,2)=1;t(2,1)=1;t(3,3)=1;
for k=1:n_proj;
    rot_matrix = t*euler_to_rot(angles(k,:));
    rot_matrices(:,:,k) = rot_matrix;
end;
% rot_matrices=rots;
% calculate inverse rotation matrices (just transpose)
inv_rot_matrices = zeros(3,3,n_proj);

%inv_rot_matrix = zeros(3);
for k=1:n_proj;
    rot_matrix = rot_matrices(:,:,k);
    inv_rot_matrix = rot_matrix'; % inv(R)=R^T
    inv_rot_matrices(:,:,k) = inv_rot_matrix;
end;