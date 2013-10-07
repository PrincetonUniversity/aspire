function [d] = q_to_d( id1, id2, q, flag )
%This is used to compute the true viewing angle difference between two
%simulated images
%   Input: 
%       id1, id2: pairs of indices to compare
%       q: quaternion
%       flag signals there is a reflection. flag==1, no reflection.
%       flag==2, reflection.
%   Output: 
%       d: relative rotations between those points
%   
%   Zhizhen Zhao Feb 10 2012
if nargin==3
    flag=1;
end;
refl=[-1, 0, 0; 0, 1, 0; 0, 0, -1];
refl=(-1)^(flag-1);
R_i=q_to_rot(q(:, id1));
R_j=q_to_rot(q(:, id2));
if flag==1
    R=R_i*R_j';
else
    R=R_i*(refl*R_j)';
end;
d=R(3,3);

end

