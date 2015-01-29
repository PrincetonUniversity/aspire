function [stack,doflip,signalmean,noisemean]=cryo_globalphaseflip(stack)
%
% CRYO_GLOBALPHASEFLIP Apply global phase flip to an image stack
%
% [stack,doflip,signalmean,noisemean]=cryo_globalphaseflip(stack)
%   Check if all images in a stack should be globally phase flipped so that
%   the molecule corresponds brigther pixels and the background corresponds
%   to darker pixels. This is done by comapring the mean in a small circle
%   around the origin (supposed to corrspond to the molecule) with the mean
%   of the noise, and making sure that the mean of the molecule is larger.
%
%   Input parameters:
%       stack  Stack of images. 
%
%   Output parameters:
%       stack       Input stack of images are flipping each image, whenever
%                   needed.
%       doflip      1 if phase was flipped and 0 otherwise.
%       signalmean  Mean of the particle part in all images.
%       noisemean   Mean of the noise in all images
%
%   Examples:
%       stack=ReadMrc('stack.mrc');
%       flippedstack=cryo_globalphaseflip(stack);
%
% Yoel Shkolnisky, January 2015

sz=size(stack);

if numel(sz)==3
    K=size(stack,3);
elseif numel(sz)==2
    K=1; % The stack has only one image.
else
    erorr('Illegal stack size');
end

if sz(1)~=sz(2)
    error('images must be suqare');
end

n=sz(1);
center=(n+1)/2;
[I,J]=meshgrid(1:n,1:n);
r=sqrt((I-center).^2+(J-center).^2);
sigind=r<round(n/4); % Indices of signal samples
noiseind=r>round(n/2*0.8);

signalmean=zeros(K,1);
noisemean=zeros(K,1);
doflip=0;

for idx=1:K
    proj=stack(:,:,idx);
    signalmean(idx)=mean(proj(sigind));
    noisemean(idx)=mean(proj(noiseind));        
end

signalmean=mean(signalmean);
noisemean=mean(noisemean);

if signalmean<noisemean
        doflip=1;
        stack=-stack;
end

