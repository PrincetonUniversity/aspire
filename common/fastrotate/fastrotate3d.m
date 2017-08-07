function OUTPUT=fastrotate3d(INPUT,R)
%FASTROTATE3D Rotate a 3D volume by a given rotation matrix.
% Input parameters:
%  INPUT    Volume to rotate, can be odd or even. 
%  R        3x3 rotation matrix.
%
% Output parameters:
%  OUTPUT   The rotated volume.
%
% Examples:
%
%   R=rand_rots(1);
%   rvol=fastrotate3d(vol,R);
%
%Yoel Shkolnisky, November 2013.

[psi,theta,phi]=rot2xyz(R);

psid=psi*180/pi;
thetad=theta*180/pi;
phid=phi*180/pi;

%figure;view3d(INPUT,0.5,'r');
tmp=fastrotate3x(INPUT,psid);
%view3d(tmp,0.5,'g');
tmp=fastrotate3y(tmp,thetad);
%view3d(tmp,0.5,'m');
OUTPUT=fastrotate3z(tmp,phid);
%view3d(OUTPUT,0.5,'b');

if isa(INPUT,'single')
    OUTPUT=single(OUTPUT);
end