function [psi,theta,phi]=rot2xyz(R)
%ROT2XYZ Convert the rotation matrix into three Euler angles.
%   Convert the rotation matrix R into three Euler angles (psi,theta,phi)
%   in x-y-z notation. That is, psi, theta, phi are rotation angles around
%   the x,y,z axes, respectively. Equivalently, R=Rz(phi)Ry(theta)Rx(psi). 
%
% Input parameters:
%  R    Rotation matrix
%
% Output parameters:
%  psi      Rotation angle around the x-axis (in radians). 
%  theta    Rotation angle around the y-axis (in radians).
%  phi      Rotation angle around the z-axis (in radians).
%
% Examples:
%
%   R=rand_rots(1);
%   [psi,theta,phi]=rot2xyz(R);
%
%Yoel Shkolnisky, November 2013.

if (1-abs(R(3,1)))>eps % R(3,1)=+/-1
	theta=-asin(R(3,1));
	psi=atan2(R(3,2)/cos(theta),R(3,3)/cos(theta));
    phi=atan2(R(2,1)/cos(theta),R(1,1)/cos(theta));
        
% % Another possible solution
%   theta=pi-theta1;
%   psi=atan2(R(3,2)/cos(theta2),R(3,3)/cos(theta2));
%   phi=atan2(R(2,1)/cos(theta2),R(1,1)/cos(theta2));
else
	phi=0; % Can be anything
	if abs(R(3,1)+1)<eps % R(3,1)=-1
		theta=pi/2;
		psi=phi+atan2(R(1,2),R(1,3));
	else
		theta=-pi/2;
		psi=-phi+atan2(-R(1,2),R(1,3));
	end
end