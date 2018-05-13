function [idx,vals_sorted]=cryo_sort_stack_bandpass(projs,lc,hc)
% CRYO_SORT_STACK_BANDPASS   Sort images by their energy in a band of frequencies.
%
% [idx,vals_sorted]=cryo_sort_stack_bandpass(projs,lc,hc)
%   Return the indices of the projections sorted by their energy in a band
%   of frequencies between lc (low cutoff) and hc (high cutoff). Both lc
%   and hc are between 0 and 0.5. Default values: lc=0.05, hc=0.2.
%
%   Returned values:
%   idx            Indices of sorted images
%   vals_sorted    Energy in the band (lc,hc) sorted from high to low.
%                  vals_sorted(i) is the energy of image idx(i).
%
% Yoel Shkolnisky, May 2018

if ~exist('lc','var')
    lc=0.05;
end

if ~exist('hc','var')
    hc=0.2;
end

fp=cfft2(projs);
L=(size(projs,1)-1)/2;
[x,y]=meshgrid(-L:L,-L:L);
idx1=find(x.^2+y.^2<(L*lc)^2);
idx2=find(x.^2+y.^2>(L*hc)^2);
vals=zeros(size(fp,3),1);
for k=1:numel(vals)
    img=fp(:,:,k);    
    img(idx2)=0;    
    if norm(img(:))>1.0e-6
        img=img-mean(img(:));
        img=img./norm(img(:)); 
        img(idx1)=0;
        v=sum(abs(img(:)).^2);
        %v=var(img(:));
    else
        v=-1;
    end
    vals(k)=v;
end
[vals_sorted,idx]=sort(vals,'descend');