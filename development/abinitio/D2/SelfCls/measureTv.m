
%Input: c1 = 180-theta, where theta is the argument for the for the
%rotation about the symmetry axis. This is since any Top View image is
%obtained by projecting the volume using a rotation about a symmetry axis
function m=measureTv(C,c1,L)
m1=measureEq2(C,c1,L);
m2=measureEq2(C,c1+90,L);
m=0.5*(m1(1)+m2(1));