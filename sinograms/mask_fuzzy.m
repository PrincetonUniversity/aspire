function [ projections,sigma ] = mask_fuzzy( projections,rad )
% Mask the projections with a radius mask of size 'rad', and compute the
% std of noise 'sigma'
% 
% Lanhui Wang, July 3, 2013

siz=size(projections,1);
center=(siz+1)/2;
m = fuzzymask(siz,2,rad,2,[center center]);
[I,J]=meshgrid(1:siz,1:siz);
r=sqrt((I-center).^2+(J-center).^2);
ind= r>rad;
n_proj=size(projections,3);
projections=reshape(projections,siz^2,n_proj);
noise=projections(ind,:);
sigma=std(noise(:));

%fuzzy mask
projections=bsxfun(@times,projections,m(:));
projections=reshape(projections,siz,siz,n_proj);
end

