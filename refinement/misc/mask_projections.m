function [projs,idx,ear]=mask_projections(rawprojs,r)
%
% Subtract the mean from all projections and maske them using a circular
% mask. 
% The Output variable 'projs' contains the masked projections. 'idx' is a
% list of the indices whose value was not masked (indices of pixels inside
% the rojecions).
%
% r is the masking radius. If r is not given, the user should enter it
% interactively based on the average of all projesions.
%
% Revisions:
%   02/03/2009 Filename changed from mask_projections_v3.m to
%              mask_projections.m
%
% Yoel Shkolnisky and Amit Singer, February 2008.

projs=zeros(size(rawprojs));
L=size(projs, 1);
if nargin<2
    mean_image=mean(rawprojs,3);
    imagesc(mean_image);
    axis image;
    r=input('masking radius?');
end

% Note: first mask and then subtract the mean, so the mean of the final
% image is zero.

[v,idx]=gaumask(size(rawprojs(:,:,1)),r);

% mask all projections
mask=zeros(size(rawprojs(:,:,1)));
mask(idx)=v;
ear=zeros(L, L, size(rawprojs, 3));
for k=1:size(rawprojs,3)
    ear(:, :, k)=(1-mask).*rawprojs(:, :, k);
    rawprojs(:,:,k) = mask.*rawprojs(:,:,k);
end
% projs=rawprojs;

%subtract mean from all projections
for k=1:size(rawprojs,3)
    p=rawprojs(:,:,k);
%    projs(:,:,k)=p-mean(p(:));
    projs(:,:,k)=p;
end


function [v,idx]=gaumask(siz,rad)
% Gaussian mask of given dimensions and radius.
center=(siz+1)./2;
[I,J]=ind2sub(siz,1:prod(siz));
r=sqrt((I-center(1)).^2+(J-center(2)).^2);
v=radialmask(r/rad);
idx=find(v>1.0e-8);
v=v(idx);
