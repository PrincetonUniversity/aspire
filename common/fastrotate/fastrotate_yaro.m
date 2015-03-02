function OUTPUT=fastrotate_yaro(INPUT,phi)

% 3-step image rotation by shearing (a copy of myrotate.m: display is disabled)
% phi is the rotation angle in grad.
% Copyright L. Yaroslavsky (yaro@eng.tau.ac.il)  
% Call OUTPUT=fastrotate(INPUT,phi); 
%
% Modified by Yoel Shkolnisky January 2011.
% Code provided by Yaroslavsky.
% Rotation is by phi degrees CW.
% Warnings were removed.


[SzX SzY] =size(INPUT);
phi=pi*phi/180;
OUTPUT=INPUT;
%subplot(2,2,1);image(OUTPUT);title('Initial image');axis image;axis off;
%drawnow;

%FIRST PASS
u=(1-cos(phi))/sin(phi+eps);
for x=1:SzX,
Ux=u*(x-SzX/2-1/2);
My=zeros(1,SzY);
r=1:SzY/2+1; Wn=1i*2*pi*Ux/SzY;
My(r)=exp(Wn*(r-1)); My(SzY:-1:SzY/2+2)=conj(My(2:SzY/2));
		spinput=fftnorm(INPUT(x,:));
		spinput=spinput.*My;
OUTPUT(x,:)=real(ifftnorm([spinput(1:SzY/2) real(spinput(SzY/2+1)) spinput(SzY/2+2:SzY)]));
end

% SECOND PASS
u=-sin(phi);
for y=1:SzY,
Uy=u*(y-SzY/2-1/2);
Mx=zeros(SzY,1);
r=1:SzX/2+1; Wn=1i*2*pi*Uy/SzX;
Mx(r)=exp(Wn*(r-1)); Mx(SzX:-1:SzX/2+2)=conj(Mx(2:SzX/2));
		spinput=fftnorm(OUTPUT(:,y));
		spinput=spinput.*Mx;
OUTPUT(:,y)=real(ifftnorm([spinput(1:SzX/2);real(spinput(SzX/2+1));spinput(SzX/2+2:SzX)]));
end

%THIRD PASS
u=(1-cos(phi))/sin(phi+eps);
for x=1:SzX,
	Ux=u*(x-SzX/2-1/2);
	My=zeros(1,SzY);
	r=1:SzY/2+1; Wn=1i*2*pi*Ux/SzY;
	My(r)=exp(Wn*(r-1)); My(SzY:-1:SzY/2+2)=conj(My(2:SzY/2));
	spinput=fftnorm(OUTPUT(x,:));
   spinput=spinput.*My;
   OUTPUT(x,:)=real(ifftnorm([spinput(1:SzY/2) real(spinput(SzY/2+1)) spinput(SzY/2+2:SzY)]));
	%OUTPUT(x,:)=fix(real(ifftnorm([spinput(1:SzY/2) real(spinput(SzY/2+1)) spinput(SzY/2+2:SzY)])));
end

%OUTPUT=(OUTPUT>=0).*OUTPUT;
%output=(OUTPUT<=255).*OUTPUT+255*ones(size(OUTPUT)).*(OUTPUT>255);