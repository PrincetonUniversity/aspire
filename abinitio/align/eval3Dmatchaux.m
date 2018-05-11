function e=eval3Dmatchaux(X,vol1,vol2)

psi=X(1);
theta=X(2);
phi=X(3);
dx=X(4);
dy=X(5);
dz=X(6);

R=Rz(phi)*Ry(theta)*Rx(psi);

vol2=fastrotate3d(vol2,R);
vol2=reshift_vol(vol2,[dx;dy;dz]);

c=corr(vol1(:),vol2(:));
c=double(c);
e=1-c;
