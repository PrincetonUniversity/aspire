function [ d, error_rot ] = check_simulation_results( class, class_refl, rot, rots_ref )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
P=size(rots_ref, 3); n_nbor=size(class, 2);
%Check class averaging results.
list=[class(:), repmat((1:P)', n_nbor, 1)];
d=zeros(P*n_nbor, 1);
rel_rot=zeros(P*n_nbor, 1);

% check acos<v_i, v_j>
for i=1:P*n_nbor
   d(i)=rot_to_d(class(i), list(i, 2), rots_ref, class_refl(i));
end;

% compute \tilde{\alpha}_{ij}
for i=1:P*n_nbor
   rel_rot(i)=rot_to_inplanerot(class(i), list(i, 2), rots_ref, class_refl(i));
end;
% compute \alpha_{ij}-\tilde{\alpha}_{ij}
error_rot=mod(rel_rot+rot(:), 360);
error_rot(error_rot>180)=error_rot(error_rot>180)-360;

end

