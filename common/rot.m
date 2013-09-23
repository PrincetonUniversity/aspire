function v=rot(u,dir,alpha)
%
% v=rot(u,dir,alpha)
%
% Rotate the vector u by an angle theta about an axis in the direction dir.
% This is implemented using Rodrigues' rotation formula.
% See http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula 
% for details.
%
% Yoel Shkolnisky, August 2010
%
ca=cos(alpha);
sa=sin(alpha);
n=cross(dir,u);
dp=dot(dir,u);
v=ca.*u+sa.*n+dp.*(1-ca).*dir;