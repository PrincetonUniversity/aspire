function cryo_normalize_background_outofcore(instackname,outstackname,r,verbose)
% CRYO_NORMALIZE_BACKGROUND_OUTOFCORE      Normalize projections
%
% cryo_normalize_background_outofcore(instackname,outstackname)
%   Normalize projections in the MRC file named instackname, and write the
%   normalized images to the MRC file outstackname. See
%   cryo_normalize_background for details of the normaliztion.
%   The function does not load the entire stack into memory.
%
% cryo_normalize_background_outofcore(instackname,outstackname,r)
%   Use pixels outside radius r as noise pixels.
%
% cryo_normalize_background_outofcore(instackname,outstackname,r,verbose)
%   Print progress messages whenever verbose is nonzero. 
%   Default: verbose=1.
%
% Example:
%   cryo_normalize_background_outofcore('instack.mrcs');
%   cryo_normalize_background_outofcore('instack.mrcs',32);
%
% Yoel Shkolnisky, May 2016.


if nargin<4
    verbose=1;
end

instack=imagestackReader(instackname);
K=instack.dim(3);
m=instack.dim(1); n=instack.dim(2);

if m~=n
    error('Images in the stack must be square.');
end

if nargin<3
    r=floor(n/2);
end

% Find indices of backgruond pixels in the images
ctr=(n+1)/2;
[I,J]=ndgrid(1:n,1:n);
radiisq=(I(:)-ctr).^2+(J(:)-ctr).^2;
background_pixels_idx=radiisq>r*r;

outstack=imagestackWriter(outstackname,instack.dim(3));

if verbose
    printProgressBarHeader;
end
for kk=1:K
    if verbose
        progressTic(kk,K);
    end
    
    proj=instack.getImage(kk);
    background_pixels=proj(background_pixels_idx);
    
    % Compute mean and standard deviation of background pixels
    mm=mean(background_pixels);
    sd=std(background_pixels);
        
    if sd<1.0e-5
        warning('Variance of background of image %d too small (sd=%5.3e). Cannot normalize...',kk,sd);
    end
    
    % Normalize the projections
    proj=(proj-mm)./sd;
    outstack.append(proj);
end
outstack.close;
