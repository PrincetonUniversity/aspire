function R=Ry(theta)
%RY Rotation matrix around the y-axis.
% Input parameters:
%  theta      Rotation angle in Radians CCW. 
%
% Output parameters:
%  R   Rotation matrix around the y-axis by an angle theta.
%
% Examples:
%  R=Ry(pi/4);
%
%Yoel Shkolnisky, November 2013.

R=[ cos(theta)  0   sin(theta);
        0       1       0;
   -sin(theta)  0   cos(theta)
    ];