function c=cryo_image_contrast_outofcore(stackname,r)
%
% CRYO_IMAGE_CONTRAST_OUTOFCORE     Estimate image contrast
%
% snr=cryo_image_contrast(stackname,r)
%   Estimate the contrast of each image in an image stack. Only pixels
%   within a radius r of the center are used to computing the contrast (to
%   eliminate noise pixels at the "corners"). The contrast of an image is
%   defined as the standard deviation of the pixels within radius r of the
%   image.
%
% snr=cryo_estimate_snr(projs)
%   Use radius r which equals to half of the image side.    
%   
% See also cryo_image_contrast
%
% Yoel Shkolnisky, December 2016.

instack=imagestackReader(stackname);
n=instack.dim(1);
if n~=instack.dim(2)
    error('Images must be square');
end

if ~exist('r','var')
    r=floor(n/2);
end

I = cart2rad(n);
idx=find(I<=r);
c=zeros(instack.dim(3),1);

for k=1:instack.dim(3)
    p=instack.getImage(k);
    c(k)=std(p(idx));
end