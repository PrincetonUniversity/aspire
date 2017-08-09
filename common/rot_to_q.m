function q=rot_to_q( rot )
%
% ROT_TO_Q 	Convert rotations to quaternions
%   
% q=rot_to_q(rot)
%   Convert a 3x3xK array of rotation matries into quaternions.
%
% See also q_to_rot
%
% Yariv Aizenbud, August 2015
% Joakim Anden, August 2017 (vectorized)

q0 = ( rot(1,1,:) + rot(2,2,:) + rot(3,3,:) + 1) / 4;
q1 = ( rot(1,1,:) - rot(2,2,:) - rot(3,3,:) + 1) / 4;
q2 = (-rot(1,1,:) + rot(2,2,:) - rot(3,3,:) + 1) / 4;
q3 = (-rot(1,1,:) - rot(2,2,:) + rot(3,3,:) + 1) / 4;


q0(q0<0) = 0;
q1(q1<0) = 0;
q2(q2<0) = 0;
q3(q3<0) = 0;

q0 = sqrt(q0);
q1 = sqrt(q1);
q2 = sqrt(q2);
q3 = sqrt(q3);

mask0 = (q0 >= q1 & q0 >= q2 & q0 >= q3);
q0(1,1,mask0) = q0(1,1,mask0) .* 1;
q1(1,1,mask0) = q1(1,1,mask0) .* sign(rot(3,2,mask0) - rot(2,3,mask0));
q2(1,1,mask0) = q2(1,1,mask0) .* sign(rot(1,3,mask0) - rot(3,1,mask0));
q3(1,1,mask0) = q3(1,1,mask0) .* sign(rot(2,1,mask0) - rot(1,2,mask0));

mask1 = (q1 >= q0 & q1 >= q2 & q1 >= q3);
q0(1,1,mask1) = q0(1,1,mask1) .* sign(rot(3,2,mask1) - rot(2,3,mask1));
q1(1,1,mask1) = q1(1,1,mask1) .* 1;
q2(1,1,mask1) = q2(1,1,mask1) .* sign(rot(2,1,mask1) + rot(1,2,mask1));
q3(1,1,mask1) = q3(1,1,mask1) .* sign(rot(1,3,mask1) + rot(3,1,mask1));

mask2 = (q2 >= q0 & q2 >= q1 & q2 >= q3);
q0(1,1,mask2) = q0(1,1,mask2) .* sign(rot(1,3,mask2) - rot(3,1,mask2));
q1(1,1,mask2) = q1(1,1,mask2) .* sign(rot(2,1,mask2) + rot(1,2,mask2));
q2(1,1,mask2) = q2(1,1,mask2) .* 1;
q3(1,1,mask2) = q3(1,1,mask2) .* sign(rot(3,2,mask2) + rot(2,3,mask2));

mask3 = (q3 >= q0 & q3 >= q1 & q3 >= q2);
q0(1,1,mask3) = q0(1,1,mask3) .* sign(rot(2,1,mask3) - rot(1,2,mask3));
q1(1,1,mask3) = q1(1,1,mask3) .* sign(rot(3,1,mask3) + rot(1,3,mask3));
q2(1,1,mask3) = q2(1,1,mask3) .* sign(rot(3,2,mask3) + rot(2,3,mask3));
q3(1,1,mask3) = q3(1,1,mask3) .* 1;

if sum(mask0|mask1|mask2|mask3) ~= size(rot, 3)
    error('Coding error.');
end

q0 = q0(:)';
q1 = q1(:)';
q2 = q2(:)';
q3 = q3(:)';

q = [q0; q1; q2; q3];

r = sqrt(sum(abs(q).^2, 1));

q = bsxfun(@times, q, 1./r);

end
