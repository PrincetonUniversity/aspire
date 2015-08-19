function [l_ij,l_ji]=commonline_R2(Ri,Rj,L)
% Another implementation of commonline_R.
% Test which is faster.
% Yoel Shkolnisky, August 2015.

Ri=Ri.';
Rj=Rj.';

Ri3=Ri(:,3);
Rj3=Rj(:,3);

%clvec=mycross(Ri3,Rj3);
clvec = [Ri3(2).*Rj3(3)-Ri3(3).*Rj3(2);
         Ri3(3).*Rj3(1)-Ri3(1).*Rj3(3);
         Ri3(1).*Rj3(2)-Ri3(2).*Rj3(1)];
     
% No need to normalize as the normalization does not affect the atan2
% below.
%
% clvec=clvec/norm(clvec);

cij=Ri.'*clvec;
cji=Rj.'*clvec;

alphaij=atan2(cij(2),cij(1));
alphaji=atan2(cji(2),cji(1));

PI=4*atan(1.0);
alphaij=alphaij+PI; % Shift from [-pi,pi] to [0,2*pi].
alphaji=alphaji+PI;

l_ij=alphaij/(2*PI)*L;
l_ji=alphaji/(2*PI)*L;

l_ij=mod(round(l_ij),L);
l_ji=mod(round(l_ji),L);
