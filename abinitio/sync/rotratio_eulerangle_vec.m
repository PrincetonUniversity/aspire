function [R,goodidx] = rotratio_eulerangle_vec(cl, k1,k2,K3s, n_theta)
% XXX FIX
% Given a 3x3 common lines matrix, where the index of each common line is
% between 1 and n_theta, compute the rotation that takes image 1 to image 2.
%
% Error codes:
%    -101 Triangle too small
%
% Yoel Shkolnisky, July 2012
% Revisions:
%   Y.S. Feb 28, 2105   Eliminate the call to ang2orth for speedup.

R=zeros(3,3,numel(K3s));
if k1==k2
    return;
end
    
TOL_idx = 1e-12;
 
idx1=cl(K3s,k2)-cl(K3s,k1); % theta3
idx2=cl(k2,K3s)-cl(k2,k1); %-theta2
idx3=cl(k1,K3s)-cl(k1,k2); % theta1

idx1=idx1(:);
idx2=idx2(:);
idx3=idx3(:);

a=cos(2*pi*idx1/n_theta); %c3
b=cos(2*pi*idx2/n_theta); %c2
c=cos(2*pi*idx3/n_theta); %c1

% Make sure that the triangle is not too small. This will happen if the
% common line between (say) cl(1,2) is close to cl(1,3). 
% To eliminate that, we require that det(G)=1+2abc-(a^2+b^2+c^2) is large
% enough.

cond=1+2.*a.*b.*c-(a.^2+b.^2+c.^2);
toosmallidx=find(cond<=1.0e-5);
goodidx=find(cond>1.0e-5);

a=a(goodidx); b=b(goodidx); c=c(goodidx);
idx2=idx2(goodidx); idx3=idx3(goodidx);
c_alpha=(a-b.*c)./sqrt(1-b.^2)./sqrt(1-c.^2);

% Fix the angles between c_ij(c_ji) and c_ik(c_jk) to be smaller than pi/2
% otherwise there will be an ambiguity between alpha and pi-alpha
ind1 = ((idx3 > n_theta/2+TOL_idx) |(idx3<-TOL_idx & idx3>-n_theta/2));
ind2 = ((idx2 > n_theta/2+TOL_idx)|(idx2<-TOL_idx & idx2>-n_theta/2) );
c_alpha((~ind1 & ind2 ) | (ind1 & ~ind2 ))=-c_alpha((~ind1 & ind2 ) | (ind1 & ~ind2 ));
    
aa=(cl(k1,k2)-1)*2*pi/n_theta;
bb=(cl(k2,k1)-1)*2*pi/n_theta;
alpha = acos(c_alpha);

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

ang1=pi-bb;
ang2=alpha;
ang3=aa-pi;
sa = sin(ang1); ca = cos(ang1);
sb = sin(ang2); cb = cos(ang2);
sc = sin(ang3); cc = cos(ang3);

R(1,1,goodidx)=cc.*ca-sc.*cb.*sa;
R(1,2,goodidx)=-cc.*sa-sc.*cb.*ca;
R(1,3,goodidx)=sc.*sb;
R(2,1,goodidx)=sc.*ca+cc.*cb.*sa;
R(2,2,goodidx)=-sa.*sc+cc.*cb.*ca;
R(2,3,goodidx)=-cc.*sb;
R(3,1,goodidx)=sb.*sa;
R(3,2,goodidx)=sb.*ca;
R(3,3,goodidx)=cb;

if ~isempty(toosmallidx)
    R(:,:,toosmallidx)=0;
end

end
