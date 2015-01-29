function R=axisangle_to_rot(rotaxis,gamma)
%
% AXISANGLE_TO_ROT  Convert axis-angle representation to rotation matrix.
%
% R=axisangle_to_rot(rotaxis,gamma)
%    Return the rotation matrix R which rotates by angle gamma (in radians)
%    around the rotation axis rotaxis. 
%
% Yoel Shkolnisky, January 2015

rotaxis=rotaxis./norm(rotaxis); % Just to be on the safe side.
N=[     0            -rotaxis(3)  rotaxis(2);
    rotaxis(3)        0         -rotaxis(1);
    -rotaxis(2)   rotaxis(1)       0];
R=expm(gamma*N);

