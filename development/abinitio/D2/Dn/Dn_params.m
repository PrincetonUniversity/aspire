
function E=Dn_params(n)
%
% Return the parameters for the Guassians that define the phantom.
% Symmetry group of the phantom: C4
%
% E has ten columns, with each column containing a different parameter for
% the ellipsoids: 
%     
%     Column 1:  A      the totla additive mass of the Gaussian
%     Column 2:  a      standard deviation in the x-direction
%     Column 3:  b      standard deviation in the y-direction
%     Column 4:  c      standard deviation in the z-direction
%     Column 5:  x0     the x-coordinate of the center of the Gaussian
%     Column 6:  y0     the y-coordinate of the center of the Gaussian
%     Column 7:  z0     the z-coordinate of the center of the Gaussian
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
%     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
%     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
%
%   The Euler angles phi, theta, and psi are given in the so-called
%   x-conventions. See http://mathworld.wolfram.com/EulerAngles.html.
%
%   For purposes of generating the phantom, the domains for the x-, y-, and 
%   z-axes span [-1,1].  Columns 2 through 7 must be specified in terms
%   of this range.
%
% Yoel Shkolnisky, June 2007
% 
% Y.S May 2013 Renamed from gaussian_params3 and comments updated.
% Eitan Rosen 2019 (Modified for Dn)

%     A      a     b     c        x0      y0      z0    phi  theta   psi
%     -------------------------------------------------------------------

%Get coordinates for regular polygon with n vertices in the xy plane
theta=2*pi/n;
thetaDeg=360/n;
vert=zeros(3,n);
down_vert1=zeros(3,n);
down_vert2=zeros(3,n);
up_vert1=zeros(3,n);
up_vert2=zeros(3,n);
dtheta=pi/18;
for j=1:n
    vert(:,j)=[cos(theta*j);sin(theta*j);0];
    down_vert1(:,j)=[cos(theta*j+dtheta);sin(theta*j+dtheta);0];
    down_vert2(:,j)=[cos(theta*j+2*dtheta);sin(theta*j+2*dtheta);0];
    up_vert1(:,j)=[cos(theta*j-dtheta);sin(theta*j-dtheta);0];
    up_vert2(:,j)=[cos(theta*j-2*dtheta);sin(theta*j-2*dtheta);0];
end
rad=0.22;
r=0.2;
%up_vert=r*vert;
%up_vert(3,:)=up_vert(3,:)+rad;
down_vert1=r*down_vert1;
down_vert1(3,:)=down_vert1(3,:)-rad;
down_vert2=r*down_vert2;
down_vert2(3,:)=down_vert2(3,:)-2*rad;
up_vert1=r*up_vert1;
up_vert1(3,:)=up_vert1(3,:)+rad;
up_vert2=r*up_vert2;
up_vert2(3,:)=up_vert2(3,:)+2*rad;
vert=r*vert;


E=zeros(5*n,10);
for j=1:n
    E(j,:)=[1,1/12,1/8,1/12,up_vert1(1,j),up_vert1(2,j),up_vert1(3,j),thetaDeg*j,0,0];
    E(n+j,:)=[1,1/12,1/8,1/12,down_vert1(1,j),down_vert1(2,j),down_vert1(3,j),thetaDeg*j,0,0];
    E(2*n+j,:)=[1,1/12,1/8,1/12,up_vert2(1,j),up_vert2(2,j),up_vert2(3,j),thetaDeg*j,0,0];
    E(3*n+j,:)=[1,1/12,1/8,1/12,down_vert2(1,j),down_vert2(2,j),down_vert2(3,j),thetaDeg*j,0,0]; 
    E(4*n+j,:)=[1,1/12,1/8,1/12,vert(1,j),vert(2,j),vert(3,j),thetaDeg*j,0,0];    
end
center=[1,1/12,1/12,1/12,0,0,0,0,0,0];
E=[E;center];
% d=0.2;
% shift_y=0.03;
% shift_z=-0.03;
% shift_x=-0.05;
% 
% shift_y2=0.2;
% shift_z2=0.4;
% shift_x2=0.02;


% E = [  
%       1     1/12  1/8 1/12      -d+shift_x    -d-shift_y      d-shift_z    0     0      0;
%       1     1/12  1/8   1/12      d-shift_x      d+shift_y      d-shift_z     0     0      0;
%       1     1/12  1/8  1/12      -d+shift_x      d+shift_y      -d+shift_z      0     0      0;
%       1     1/12  1/8  1/12      d-shift_x    -d-shift_y      -d+shift_z      0     0      0;
%       
% %       1     1/12  1/12 1/12      d+shift_x2     d+shift_y2      -d+shift_z2    0     0      0;
% %       1     1/12  1/12   1/12       -d-shift_x2      -d-shift_y2      -d+shift_z2     0     0      0;
% %       1     1/12  1/12  1/12      d+shift_x2      -d-shift_y2      d-shift_z2      0     0      0;
% %       1     1/12  1/12   1/12       -d-shift_x2    d+shift_y2      d-shift_z2      0     0      0;
%      ];


