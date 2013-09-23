function R = rotratio_eulerangle(cl, n_theta)
% 
% Given a 3x3 common lines matrix, where the index of each common line is
% between 1 and n_theta, compute the rotation that takes image 1 to image 2.
%
% Error codes:
%    -101 Triangle too small
%
% Yoel Shkolnisky, July 2012

TOL_idx = 1e-12;
 
idx1=cl(3,2)-cl(3,1); % theta3
idx2=cl(2,3)-cl(2,1); %-theta2
idx3=cl(1,3)-cl(1,2); % theta1

a=cos(2*pi*idx1/n_theta); %c3
b=cos(2*pi*idx2/n_theta); %c2
c=cos(2*pi*idx3/n_theta); %c1

% Make sure that the triangle is not too small. This will happen if the
% common line between (say) cl(1,2) is close to cl(1,3). 
% To eliminate that, we require that det(G)=1+2abc-(a^2+b^2+c^2) is large
% enough.

if 1+2*a*b*c-(a^2+b^2+c^2)<1.0e-5
    R=-101;
    return;
end

c_alpha=(a-b.*c)./sqrt(1-b.^2)./sqrt(1-c.^2);

    % fix the angles between c_ij(c_ji) and c_ik(c_jk) to be smaller than pi/2
    % otherwise there will be an ambiguity between alpha and pi-alpha
ind1 = ((idx3 > n_theta/2+TOL_idx) |(idx3<-TOL_idx & idx3>-n_theta/2));
ind2 = ((idx2 > n_theta/2+TOL_idx)|(idx2<-TOL_idx & idx2>-n_theta/2) );
    
if (~ind1 && ind2 ) || (ind1 && ~ind2 )
    c_alpha = -c_alpha;
end
% c_alpha = - (cos3-cos2*cos1)/(sin2*sin1)
    
aa=(cl(1,2)-1)*2*pi/n_theta;
bb=(cl(2,1)-1)*2*pi/n_theta;
alpha = acos(c_alpha);
R =ang2orth(pi-bb, alpha, aa-pi); % the relative rotation can be represented using ZYZ Euler angle (pi-bb, alpha, aa-pi).
end

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

   