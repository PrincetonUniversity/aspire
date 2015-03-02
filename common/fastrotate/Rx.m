function R=Rx(psi)
%RX Rotation matrix around the x-axis.
% Input parameters:
%  psi      Rotation angle in Radians CCW. 
%
% Output parameters:
%  R   Rotation matrix around the x-axis by an angle psi.
%
% Examples:
%  R=Rx(pi/4);
%
%Yoel Shkolnisky, November 2013.

R=[ 1       0           0;
    0   cos(psi)    -sin(psi);
    0   sin(psi)     cos(psi)
    ];