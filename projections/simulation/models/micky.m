function E=micky
% Return the parameters for the Guassians that define the phantom.
% Symmetry group of the phantom: C1 (no symmetry).
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
% Yoel Shkolnisky August 2013
% Phantom for molecule reconstruction experiments (micky)
E =    [ .3  .500   .500  .500    0     0     0       0      0      0
         .3  .050   .300  .300    0    .40   .40      0      0      0   %ears
         .3  .050   .300  .300    0   -.40   .40      0      0      0   %ears
         .3  .300   .100  .100   .50    0     0       0      0      0
         .3  .120   .120  .150   .37  -.20   .27      0      0      0
         .3  .120   .120  .150   .37   .20   .27      0      0      0
         .3  .070   .200  .200   .50    0    .10     90     45    -90   %mouth
         .3  .100   .100  .100    0    .70     0      0      0      0   %symmetry point
        ];
