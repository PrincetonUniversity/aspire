
%Input:  ax = unit vector in R^3
%        t = angle in radians
%Output: R=rotation by ang about axis
function R=axisang2rot(ax,t)
K=[0,-ax(3),ax(2);
   ax(3),0,-ax(1);
   -ax(2),ax(1),0];
R=eye(3)+sin(t)*K+(1-cos(t))*K*K;
    
