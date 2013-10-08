function M=fastrotateprecomp(SzX,SzY,phi)
%
% Compute the interpolation tables required to rotate an image with SzX
% rows and SzY columns by an angle phi CCW.
%
% This function is used to accelerate fastrotate, in case many images are
% needed to be rotated by the same angle. In such a case it allows to
% precompute the interpolation tables only once instead of computing them
% for each image.
%
% M is a structure containing phi, Mx, My, where Mx and My are the
% interpolation tables used by fastrotate.
%
% Yoel Shkolnisky, January 2011.

% Adjust the rotation angle to be between -45 and 45 degrees.
[phi,mult90]=adjustrotate(phi);

phi=pi*phi/180;
phi=-phi; % To match Yaroslavsky's code which rotates CW.

if mod(SzY,2)==0
    cy=SzY/2+1;
    sy=1/2; % By how much should we shift the cy to get the center of the image
else
    cy=(SzY+1)/2;
    sy=0;
end

if mod(SzX,2)==0
    cx=SzX/2+1; % By how much should we shift the cy to get the center of the image
    sx=1/2;
else
    cx=(SzX+1)/2;
    sx=0;
end

% Precompte My and Mx
My=zeros(SzY,SzX);
r=1:cy;
u=(1-cos(phi))/sin(phi+eps);
alpha1=2*pi*1i*(r-1)./SzY;
for x=1:SzX
    Ux=u*(x-cx+sx);
    My(r,x)=exp(alpha1.*Ux);
    My(SzY:-1:cy+1,x)=conj(My(2:cy-2*sy,x));
end
My=My.'; % Remove when implementing using the loops below.

Mx=zeros(SzX,SzY);
r=1:cx;
u=-sin(phi);
alpha2=2*pi*1i*(r-1)./SzX;
for y=1:SzY,
    Uy=u*(y-cy+sy);
    Mx(r,y)=exp(alpha2.*Uy); Mx(SzX:-1:cx+1,y)=conj(Mx(2:cx-2*sx,y));
end

M.phi=phi;
M.Mx=Mx;
M.My=My;
M.mult90=mult90;
