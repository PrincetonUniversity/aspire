function R=Rz(phi)
%RZ Rotation matrix around the z-axis.
% Input parameters:
%  phi      Rotation angle in Radians CCW. 
%
% Output parameters:
%  R   Rotation matrix around the z-axis by an angle phi.
%
% Examples:
%  R=Rz(pi/4);
%
%Yoel Shkolnisky, November 2013.

R=[ cos(phi)  -sin(phi) 0;
    sin(phi)   cos(phi) 0; 
        0          0    1
    ];