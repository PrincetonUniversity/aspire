function OUTPUT=fastrotate_ref(INPUT,phi)
%
% 3-step image rotation by shearing.
% Rotate INPUT by phi degrees CCW.
% 
% Input parameters:
%  INPUT    Image to rotate, can be odd or eve.
%  phi      Rotation angle in degrees.
%
% Output parameters:
%  OUTPUT   The rotated image.
%
% Based on code by L. Yaroslavsky (yaro@eng.tau.ac.il). 
% The following changes have been implemented
%   1. Adapt to both odd and even sizes.
%   2. Avoid calling to fftnorm and ifftnorm. The overhead is too large.
%
% Yoel Shkolnisky, January 2011


[SzX SzY] =size(INPUT);
phi=-pi*phi/180; % The minus is since Yaroslavsky's original code rotates CW.
OUTPUT=INPUT;

% % Determine centers for even and odd sizes
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
% 
%INPUT=INPUT./SzY;
%sqtN=sqrt(SzY);


%FIRST PASS
u=(1-cos(phi))/sin(phi+eps);
My=zeros(1,SzY);
r=1:cy;
alpha1=2*pi*1i*(r-1)./SzY;
for x=1:SzX,
    Ux=u*(x-cx+sx);    
%    Wn=1i*2*pi*Ux/SzY;
%    My(r)=exp(Wn*(r-1)); My(SzY:-1:SzY/2+2)=conj(My(2:SzY/2));
%    My(r)=exp(alpha1.*Ux); My(SzY:-1:halfSzY+2)=conj(My(2:halfSzY));
    My(r)=exp(alpha1.*Ux); My(SzY:-1:cy+1)=conj(My(2:cy-2*sy));
    spinput=fft(INPUT(x,:));
    spinput=spinput.*My;
    %    OUTPUT(x,:)=real(ifftnorm([spinput(1:cy-1) real(spinput(cy)) spinput(cy+1:SzY)]));
    %    OUTPUT(x,:)=real(ifft([spinput(1:cy-1) real(spinput(cy)) spinput(cy+1:SzY)]));
    OUTPUT(x,:)=real(ifft(spinput)); % The above real() sees redundant.
    %    OUTPUT(x,:)=OUTPUT(x,:).*sqtN;
end

% SECOND PASS
u=-sin(phi);
Mx=zeros(SzY,1);
r=1:cx;
alpha2=2*pi*1i*(r-1)./SzX;
for y=1:SzY,
    Uy=u*(y-cy+sy);    
%    Wn=1i*2*pi*Uy/SzX;
%    Mx(r)=exp(Wn*(r-1)); Mx(SzX:-1:SzX/2+2)=conj(Mx(2:SzX/2));
%    Mx(r)=exp(alpha2.*Uy); Mx(SzX:-1:halfSzX+2)=conj(Mx(2:halfSzX));
    Mx(r)=exp(alpha2.*Uy); Mx(SzX:-1:cx+1)=conj(Mx(2:cx-2*sx));
    spinput=fft(OUTPUT(:,y));
    spinput=spinput.*Mx;
    %    OUTPUT(:,y)=real(ifftnorm([spinput(1:cx-1);real(spinput(cx));spinput(cx+1:SzX)]));
    %     OUTPUT(:,y)=real(ifft([spinput(1:cx-1);real(spinput(cx));spinput(cx+1:SzX)]));
    OUTPUT(:,y)=real(ifft(spinput)); % The above real() sees redundant.
end

%THIRD PASS
u=(1-cos(phi))/sin(phi+eps);
My=zeros(1,SzY);
r=1:cy;
for x=1:SzX,
    Ux=u*(x-cx+sx);    
%    Wn=1i*2*pi*Ux/SzY;
%    My(r)=exp(Wn*(r-1)); My(SzY:-1:SzY/2+2)=conj(My(2:SzY/2));
%    My(r)=exp(alpha1.*Ux); My(SzY:-1:halfSzY+2)=conj(My(2:halfSzY));
    My(r)=exp(alpha1.*Ux); My(SzY:-1:cy+1)=conj(My(2:cy-2*sy));
    spinput=fft(OUTPUT(x,:));
    spinput=spinput.*My;
    %    OUTPUT(x,:)=real(ifftnorm([spinput(1:cy-1) real(spinput(cy)) spinput(cy+1:SzY)]));
    %    OUTPUT(x,:)=OUTPUT(x,:).*sqtN;
    %OUTPUT(x,:)=fix(real(ifftnorm([spinput(1:SzY/2) real(spinput(SzY/2+1)) spinput(SzY/2+2:SzY)])));
    %OUTPUT(x,:)=real(ifft([spinput(1:cy-1) real(spinput(cy)) spinput(cy+1:SzY)]));
    OUTPUT(x,:)=real(ifft(spinput)); % The above real() sees redundant.
end

%OUTPUT=(OUTPUT>=0).*OUTPUT;
%output=(OUTPUT<=255).*OUTPUT+255*ones(size(OUTPUT)).*(OUTPUT>255);