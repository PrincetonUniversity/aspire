function [rel_rot] = rot_to_inplanerot( id1, id2, rots, flag )
%Estimate the best in-plane rotation from the 3x3 matrices. This is used to
%compute the true relative in-plane rotation.
%   Input: 
%       id1, id2: pairs of indices to compare
%       q: quaternion
%       flag: flag==1, no reflection. flag==2, reflection.
%   Output: 
%       rel_rot: relative rotations between those points
%   
%   Zhizhen Zhao Feb 10 2012

if nargin==3;
    flag=1;
end;
refl=[1, 0, 0; 0, -1, 0; 0, 0, -1];
R_i=rots(:, id1);
R_j=rots(:, id2);
if flag==1
    R=R_i*R_j';
else
    R=(refl*R_i)*R_j';
end;
rel_rot=atan2(R(2,1)-R(1,2), R(1,1)+R(2,2)); %Use the formula in Viewing angle classification 
rel_rot=rel_rot*180/pi;

end

