function q=rot_to_q( rot )
%
% ROT_TO_Q 	Convert rotaion to quaternion
%   
% q=rot_to_q(rot)
%   Convert a single 3x3 rotation matrix into a quaternion.
%
% See also q_to_rot
%
% Yariv Aizenbud, August 2015


q0 = ( rot(1,1) + rot(2,2) + rot(3,3) + 1) / 4;
q1 = ( rot(1,1) - rot(2,2) - rot(3,3) + 1) / 4;
q2 = (-rot(1,1) + rot(2,2) - rot(3,3) + 1) / 4;
q3 = (-rot(1,1) - rot(2,2) + rot(3,3) + 1) / 4;
if(q0 < 0)
    q0 = 0;
end
if(q1 < 0) 
    q1 = 0;
end
if(q2 < 0) 
    q2 = 0;
end

if(q3 < 0) 
    q3 = 0;
end
q0 = sqrt(q0);
q1 = sqrt(q1);
q2 = sqrt(q2);
q3 = sqrt(q3);
if(q0 >= q1 && q0 >= q2 && q0 >= q3) 
    q0 = q0 * 1;
    q1 = q1 * sign(rot(3,2) - rot(2,3));
    q2 = q2 * sign(rot(1,3) - rot(3,1));
    q3 = q3 * sign(rot(2,1) - rot(1,2));
elseif(q1 >= q0 && q1 >= q2 && q1 >= q3) 
    q0 = q0 * sign(rot(3,2) - rot(2,3));
    q1 = q1 * 1;
    q2 = q2 * sign(rot(2,1) + rot(1,2));
    q3 = q3 * sign(rot(1,3) + rot(3,1));
elseif(q2 >= q0 && q2 >= q1 && q2 >= q3)
    q0 = q0 * sign(rot(1,3) - rot(3,1));
    q1 = q1 * sign(rot(2,1) + rot(1,2));
    q2 = q2 * +1;
    q3 = q3 * sign(rot(3,2) + rot(2,3));
elseif(q3 >= q0 && q3 >= q1 && q3 >= q2)
    q0 = q0 * sign(rot(2,1) - rot(1,2));
    q1 = q1 * sign(rot(3,1) + rot(1,3));
    q2 = q2 * sign(rot(3,2) + rot(2,3));
    q3 = q3 * +1;
else
    printf('coding error\n');
end
r = norm([q0, q1, q2, q3]);
q0 = q0/r;
q1 = q1/r;
q2 = q2/r;
q3 = q3/r;

q = [q0;q1;q2;q3];
end