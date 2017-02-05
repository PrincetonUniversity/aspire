function c=cryo_image_contrast(projs,r)
%
% CRYO_IMAGE_CONTRAST     Estimate image contrast
%
% snr=cryo_image_contrast(projs,r)
%   Estimate the contrast of each image in an image stack. Only pixels
%   within a radius r of the center are used to computing the contrast (to
%   eliminate noise pixels at the "corners"). The contrast of an image is
%   defined as the standard deviation of the pixels within radius r of the
%   image.
%
% snr=cryo_image_contrast(projs)
%   Use radius r which equals to half of the image side.    
%   
% Example: 
%   % Show 5 images with highest contrast and 5 images with the lowest
%   % contrast
%   projs=ReadMRC('projs.mrc');
%   c=cryo_image_contrast(projs);
%   [~,ii]=sort(c);
%   viewstack(projs(:,:,[ii(1:5),ii(end-4:end)]),5,2)
%
% Yoel Shkolnisky, December 2016.


if ~isnumeric(projs)
    error('First argument must be a 2 or 3 dimensional array');
end

n=size(projs,1);
if n~=size(projs,2)
    error('Images must be square');
end

if ~exist('r','var')
    r=floor(n/2);
end

I = cart2rad(n);
idx=find(I<=r);
c=zeros(size(projs,3),1);

for k=1:size(projs,3)
    p=projs(:,:,k);
    c(k)=std(p(idx));
end
