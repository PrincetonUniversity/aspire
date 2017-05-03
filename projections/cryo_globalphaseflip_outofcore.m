function [doflip,signalmean,noisemean]=cryo_globalphaseflip_outofcore(instackname,outstackname)
%
% CRYO_GLOBALPHASEFLIP_OUTOFCORE Apply global phase flip to an image stack
%
% [doflip,signalmean,noisemean]=cryo_globalphaseflip_outofcore(instackname,outstackname)
%   Check if all images in a stack should be globally phase flipped so that
%   the molecule corresponds brigther pixels and the background corresponds
%   to darker pixels. This is done by comapring the mean in a small circle
%   around the origin (supposed to corrspond to the molecule) with the mean
%   of the noise, and making sure that the mean of the molecule is larger.
%   Input images are read from the MRC file instackname and are written in
%   MRC file outstackname.
%
% [doflip,signalmean,noisemean]=cryo_globalphaseflip_outofcore(instackname)
%   No output images are written to disk, but only output parameters are
%   computed.
%
%   The function applies global phaseflip without loading all images into
%   memory.
%
%   Input parameters:
%       instackname   Name of MRC file containing input images.
%       outstackname  Name of MRC file into which corrected image are
%                     written.
%
%   Output parameters:
%       doflip      1 if phase was flipped and 0 otherwise.
%       signalmean  Mean of the particle part in all images.
%       noisemean   Mean of the noise in all images
%
%   Examples:     
%       cryo_globalphaseflip_outofcore('instack.mrc','outstack.mrc');
%
%   See also: cryo_globalphaseflip
%
% Yoel Shkolnisky, May 2016

instack=imagestackReader(instackname,100);
sz=instack.dim;

if numel(sz)==3
    K=instack.dim(3);
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
    proj=instack.getImage(idx);
    signalmean(idx)=mean(proj(sigind));
    noisemean(idx)=mean(proj(noiseind));        
end

signalmean=mean(signalmean);
noisemean=mean(noisemean);

if signalmean<noisemean
        doflip=1;        
end

if nargin==2
    % Also write the flipped projections to the output file
    outstack=imagestackWriter(outstackname,K,1,100);
    flag=(-2)*doflip+1; % Sign multiplier for the images
    for idx=1:K
        proj=instack.getImage(idx);
        proj=proj.*flag;
        outstack.append(proj);
    end
    outstack.close;    
end