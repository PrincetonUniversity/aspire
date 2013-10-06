function mr=grotate(m,theta,kernelsize)
% function [mr,rft]=grotate(m,theta[,kernelsize]);
% 2D rotation using the fast gridding algorithm of Penczek.
% The image is padded by 1.25x and the interpolation is done in Fourier
% space.
% m is the input 2d image, n x n with n divisible by 4.  mr is the rotated
% output image.
% m is assumed to be strictly band-limited, i.e. containing no spatial
% frequencies above 0.5-1/n in any direction.
%
% Kernelsize is an optional argument.  It is an odd number: 3, 5, 7 or 9.
% The default value is 3.

% F. Sigworth 16 May 07
% Vectorized code.
% takes 100 ms to rotate a 48x48 image.  Half this time is in the
% interpolation, the rest is in the various FTs and other operations.


% default argument
if nargin<3
    kernelsize=3;
end;
nw=kernelsize;

% if no rotation is to be done
if theta==0
    mr=m;
end;

[n n1]=size(m);

[w1 ov]=gridMakeKaiserTable(kernelsize,'grid');
w1=w1';

wcomp=gridMakePreComp(n,nw);
P=gridMakePaddedFT(m,'grid',wcomp);
P2=gridMakeNullFT(n,2);

% Interpolation is done here.
np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;
sp1=P.sp1;
ovctr=ov/2+1;

% Vectorized code.
plane=zeros(np*np,1);  % vector to receive the data
nw2=floor(kernelsize/2);  % half-width of kernel (=(nw-1)/2).
s=sin(theta);
c=cos(theta);

% Source coordinates
rmax=np/2+1;
[is,js]=ndgrid(-np/2:np/2-1);
is=reshape(is,np*np,1);
js=reshape(js,np*np,1);

% Object coordinates
ip=c*is+s*js+np1ctr;
ip=min(sp1+np,max(ip,sp1+1));  % prevent coords from going out of bounds
ipint=round(ip);
ipfrac=floor((ip-ipint)*ov)+ovctr;

jp=-s*is+c*js+np1ctr;
jp=min(sp1+np,max(jp,sp1+1));
jpint=round(jp);
jpfrac=floor((jp-jpint)*ov)+ovctr;

addrs=ipint+np1*(jpint-1);  % 1-dim addresses in the PadFT array.
for i1=-nw2:nw2
    for j1=-nw2:nw2
        plane=plane+w1(ipfrac,i1+nw2+1).*w1(jpfrac,j1+nw2+1)...
            .*P.PadFT(addrs+i1+(np1*j1));
    end;
end;
% % Convert the vector to a 2d array
P2.PadFT(sp1+1:sp1+np,sp1+1:sp1+np)=reshape(plane,np,np);
mr=gridRecoverRealImage(P2);
