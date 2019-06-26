
function E=D2_params
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


%     A      a     b     c        x0      y0      z0    phi  theta   psi
%     -------------------------------------------------------------------
d=0.12;
shift_y=0.03;
shift_z=-0.03;
shift_x=-0.05;

shift_y2=0.2;
shift_z2=0.4;
shift_x2=0.02;


E = [  
      1     1/12  1/8 1/12      -d+shift_x    -d-shift_y      d-shift_z    0     0      0;
      1     1/12  1/8   1/12      d-shift_x      d+shift_y      d-shift_z     0     0      0;
      1     1/12  1/8  1/12      -d+shift_x      d+shift_y      -d+shift_z      0     0      0;
      1     1/12  1/8  1/12      d-shift_x    -d-shift_y      -d+shift_z      0     0      0;
      
%       1     1/12  1/12 1/12      d+shift_x2     d+shift_y2      -d+shift_z2    0     0      0;
%       1     1/12  1/12   1/12       -d-shift_x2      -d-shift_y2      -d+shift_z2     0     0      0;
%       1     1/12  1/12  1/12      d+shift_x2      -d-shift_y2      d-shift_z2      0     0      0;
%       1     1/12  1/12   1/12       -d-shift_x2    d+shift_y2      d-shift_z2      0     0      0;
     ];


