function [stack,mean_bg,sd_bg]=cryo_normalize_background(stack,r,verbose)
%
% CRYO_NORMALIZE_BACKGROUND Normalize background to mean 0 and std 1.
%
% stack=cryo_normalize_background(stack,r)
%   Estimate the mean and std of each image in the stack using pixels
%   outside radius r (pixels), and normalize the image such that the
%   background has mean 0 and std 1. Each image in the stack is corrected
%   separately.
%
% stack=cryo_normalize_background(stack)
%   Set r to be half the image side.
%
% Example:
% stack2=cryo_normalize_background(stack,55);
%
% Yoel Shkolnisky, February 2015
%
% Revisions:
% Y.S.  August 2015     Add default radis and verbose flag.

if nargin<3
    verbose=1;
end

K=size(stack,3);
m=size(stack,1); n=size(stack,2);

if m~=n
    error('Images in the stack must be square.');
end

if nargin<2
    r=floor(n/2);
end

% Find indices of backgruond pixels in the images
ctr=(n+1)/2;
[I,J]=ndgrid(1:n,1:n);
radiisq=(I(:)-ctr).^2+(J(:)-ctr).^2;
background_pixels_idx=radiisq>r*r;

if verbose
    printProgressBarHeader;
end
sd_bg = zeros(1,K);
mean_bg = zeros(1,K);
for kk=1:K
    if verbose
        progressTic(kk,K);
    end
    
    proj=stack(:,:,kk);
    background_pixels=proj(background_pixels_idx);
    
    % Compute mean and standard deviation of background pixels
    mm=mean(background_pixels);
    sd=std(background_pixels);
%    sprintf('Subtracting background mean: %f and normalizing by background stdev: %f',mm, sd)
    % Normalize the projections
    
    if sd<1.0e-5
        warning('Variance of background of image %d too small (sd=%5.3e). Cannot normalize...',kk,sd);
    end
    
    proj=(proj-mm)./sd;
    stack(:,:,kk)=proj;
    
    sd_bg(kk) = sd;
    mean_bg(kk) = mm;
end
