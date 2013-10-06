function [mr,rft]=grotate1(m,theta,kernelsize)
% function [mr,rft]=grotate(m,theta,kernelsize);
% 2D rotation using the fast gridding algorithm of Penczek.
% The image is padded by 1.25x and the interpolation is done in Fourier
% space.
% m is the input 2d image, n x n with n divisible by 4.  rm is the rotated
% output image.
% m is assumed to be strictly band-limited, i.e. containing no spatial
% frequencies above 0.5-1/n in any direction.
%
% Kernelsize is an optional argument.  It is an odd number: 3, 5, 7 etc.
% The default value is 3.  It must be smaller than n/10.
% (rft is a debugging variable.  It is the rotated fourier transform of the
% object.)

% F. Sigworth 26 July 06
% Sped up the interpolation by using vector product fs 28 July 06
% Added padding for the ft 4 Aug 06
% calls: fuzzymask, kaiser
% takes ~0.4 s for a 64x64 image.

% parameters and static variables
AlphaValues=[0 0 5.2 0 10.2 0 13 0 17];  % nice alpha values for the oversampling we use.
ov=1024;  % oversampling of window for look-up
persistent w1 % Cache for the w1 lookup table.

% default argument
if nargin<3
    kernelsize=3;
end;

nw=kernelsize;
alpha=gridGetAlpha(nw);
if alpha==0
    error('Invalid kernel size; it must be odd.');
end;


if (numel(w1)~=nw*ov) % if the cached table is absent or not the right size
    % Create the interpolation kernel
    dw=1/ov;  % fractional increment.  We assume ov is even
    k=(-nw/2+dw/2:dw:nw/2-dw/2)';  % nw*ov space points
    w=kaiser(nw/2,alpha,k);  % Compute the Kaiser window
    w=w*ov/sum(w);  % normalize it.
    w1=zeros(nw,ov);

    % Make the 1D lookup table w1(:,i).  i=1 corresponds to a shift between -0.5 and
    % 0.5+dw, so we assign it the mean shift -0.5+dw/2.
    % i=ov corresponds to a mean shift of 0.5-dw/2.
    for i=1:ov
        w1(:,i)=w(ov-i+1:ov:nw*ov);
    end;
end;


[n n1]=size(m);

wcomp=gridMakePreComp(n,nw);

P=gridMakePaddedFT(m,'grid',wcomp);

np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;

sp=P.sp;
sp1=P.sp1;  % shift of origin going from np to np1
np1ctr=n/2+1+sp+sp1;


  % Create the interpolation kernel
if (numel(w1)~=nw*ov) % if the cached table is absent or not the right size
  dw=1/ov;  % fractional increment.  We assume ov is even
  k=(-nw/2+dw/2:dw:nw/2-dw/2)';  % nw*ov space points
  w=kaiser(nw/2,alpha,k);  % Compute the Kaiser window
  w=w*ov/sum(w);  % normalize it.
  w1=zeros(nw,ov);

  % Make the 1D lookup table w1(:,i).  i=1 corresponds to a shift between -0.5 and
  % 0.5+dw, so we assign it the mean shift -0.5+dw/2.
  % i=ov corresponds to a mean shift of 0.5-dw/2.
  for i=1:ov
    w1(:,i)=w(ov-i+1:ov:nw*ov);
  end;
end;


% Interpolation is done here.

s=sin(theta);
c=cos(theta);
rot=[c s;-s c];  % rotation matrix

r=round(np/2+1);  % radius of fourier-space region to consider.
% We extend it beyond the original n so that we can roll off the Fourier space image.

centerp=[np1ctr np1ctr]';  % coordinate of the origin
rft1=zeros(np1,np1); % array to receive the padded, rotated image
nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).

for i=np1ctr-r:np1ctr+r;  % Everything beyond r is assumed to be zero:
  %   we compute points only in the disc of radius r
  jextent=ceil(sqrt(r.^2-(i-np1ctr).^2));
  for j=np1ctr-jextent:np1ctr+jextent  % loop over each output point!
    p=rot*([i j]'- centerp)+centerp;
    pint=round(p);
    pfrac=floor((p-pint)*ov)+ov/2+1;
    i1=pint(1);
    j1=pint(2);
    rft1(i,j)=w1(:,pfrac(1))'*P.PadFT(i1-nw2:i1+nw2,j1-nw2:j1+nw2)*w1(:,pfrac(2));
  end;
end;

% % Remove padding and transform back to real space.
% rft=rft1(sp1+1:sp1+np,sp1+1:sp1+np);
% rm=fftshift(real(ifftn(fftshift(rft))));
% 
% rmcm=rm.*disc(np,n/2+1,[np/2+1,np/2+1]);  % mask out the outlying part
% mr=rmcm(1+sp:n+sp,1+sp:n+sp);  % extract the un-padded result.
P.PadFT=rft1;
mr=gridRecoverRealImage(P);
