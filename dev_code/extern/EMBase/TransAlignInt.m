function [aligned,shifts,amps]=TransAlignInt(stack,ref,ccprior,maxshift)
% function [aligned,shifts]=TransAlignInt(stack,ref[,ccmask,maxshift]) 
% Perform a translational alignment of each image in the stack to the
% reference. The stack can also be a single image. Only integer shifts are
% performed, and they are performed circularly. This routine is used to
% make a rough alignment of raw images. The optional ccmask argument is
% a realspace mask that is added to the cc function.  This is typically
% a quadratic to restrict the centering range according to ML principles.
% According to Sigworth (98) it should be equal to sigma^2 x ln( f(trans))
% where sigma is the noise sd and f is the pdf of translations.  For a
% Gaussian translation error of sd sigma_trans pixels, 
% ccprior = -sigma^2 R^2/(2 sigma_trans^2) where R=Radius(n). 
% The default maxshift is n/4 pixels.
% The returned vector amps are the values of the cc peaks, and shifts is a
% 2 x nim array of integer values, the shifts applied to make the aligned 
% images.  alignedImage=circshift(image,shift);

if (ndims(stack)<3)
    ni=1;
    [n1 n2]=size(stack);
else
    [n1 n2 ni]=size(stack);
end;
if nargin<4
    maxshift=n1/4;
end;
if nargin<3
    ccprior=0;
end;

fref=conj(fftn(ref));
aligned=stack;
shifts=zeros(2,ni);
amps=zeros(ni,1);

for i=1:ni
    x=stack(:,:,i);
    fx=fftn(x);
    cc=fftshift(real(ifftn(fref.*fx)));
    [mx dx dy]=max2d(cc+ccprior);
    shift=[n1/2+1-dx n2/2+1-dy];
    shift=max(-maxshift,min(maxshift,shift));  % Not needed due to mask.
    aligned(:,:,i)=circshift(x,shift);
    shifts(:,i)=shift';
    amps(i)=mx;
end;
