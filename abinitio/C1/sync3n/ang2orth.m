
function orthm = ang2orth(ang1,ang2,ang3)
% Convert the Euler angles with ZXZ conversion to rotation matrices
% Euler angle (a,b,c) to rotation
% ra = [  ca,  -sa,  0; ...
%     sa,  ca,  0; ...
%     0,   0,  1];
% rb = [  1,   0,   0; ...
%     0,   cb, -sb;...
%     0,  sb, cb];
% rc = [  cc,  -sc,  0; ...
%     sc,  cc,  0; ...
%     0,   0,  1];
% orthm = rc*rb*ra;
% ca is short for cos(a) and sa is for sin(a).
%
% This function does the conversion simultanously for N Euler angles.
%
% Input:
%      ang1/ang2/ang3: a vector of length N. (ang1(i),ang2(i),ang3(i))
%      determines a rotation matrix orthm(:,:,i)
% Ouput:
%      orthm: rotation matrixs of size 3x3xN
%
ang1=ang1(:);
ang2=ang2(:);
ang3=ang3(:);
sa = sin(ang1); ca = cos(ang1);
sb = sin(ang2); cb = cos(ang2);
sc = sin(ang3); cc = cos(ang3);
n=length(ang1);
orthm=zeros(3,3,n);
orthm(1,1,:)=cc.*ca-sc.*cb.*sa;
orthm(1,2,:)=-cc.*sa-sc.*cb.*ca;
orthm(1,3,:)=sc.*sb;
orthm(2,1,:)=sc.*ca+cc.*cb.*sa;
orthm(2,2,:)=-sa.*sc+cc.*cb.*ca;
orthm(2,3,:)=-cc.*sb;
orthm(3,1,:)=sb.*sa;
orthm(3,2,:)=sb.*ca;
orthm(3,3,:)=cb;
end
