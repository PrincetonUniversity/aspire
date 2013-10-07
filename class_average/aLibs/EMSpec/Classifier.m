function [means nmembers indices fmeans eigim svals]=Classifier(img,nclasses,nfactors,decim,cmask)
% function [means nmembers indices fmeans eigim svals]=Classifier(img,nclasses,nfactors,decim,cmask)
% Take the set of masked images and perform SVD using nfactors dimensions,
% and then kmeans clustering to yield nclasses.
% nclasses=100;
% nfactors=40;
% decim=2;
% if nargin<4
%     decim=2;
% end;
maskThreshold=0.1;

[nx0 ny0 nimgs]=size(img);

mimgs=img;  % allocate it.
% Make masked images
if nargin >4  % we have a mask
    for i=1:nimgs
        mimgs(:,:,i)=img(:,:,i).*cmask;
    end;
else
    cmask=1;
end;
% Make downsampled images for svd purposes
if decim~=1
    nx=nx0/decim;
    ny=ny0/decim;
    % Fourier downsampling
    mimgs=Downsample(mimgs,nx,1);  % I assume square images here.
    cmask=Downsample(cmask,nx);
else
    nx=nx0;
    ny=ny0;
end;
pix=cmask>maskThreshold;  % We will extract just these pixels.
imgs1d=zeros(sum(pix(:)),nimgs);  % Images as 1D vectors
for i=1:nimgs
    im=mimgs(:,:,i);
    imgs1d(:,i)=im(pix);  % each image is a column
end;

disp('svd...');
tic
[u s v]=svd(imgs1d);
% [u s v]=svds(imgs1d,nfactors);
disp('done.');
toc
% u=u(:,1:nfactors);
% s=s(1:nfactors,1:nfactors);
% v=v(:,1:nfactors);
% u has eigenimages as columns, nx*ny by ndims
% s*v' is nfactors by nimgs, the vectors.
factors=s*v';
eigim=zeros(nx,ny,nfactors);
eig1=zeros(nx,ny);
for i=1:nfactors  % Restore the pixels to images.
    eig1(pix)=u(:,i);
    eigim(:,:,i)=eig1;
end;
% eigim=reshape(u, nx, ny, nfactors);
svals=diag(s);
truncFactors=factors(1:nfactors,:);
disp('kmeans...');
tic
[indices fmeans]=kmeans(truncFactors',nclasses);
disp('done.');
toc

% We could get low-resolution means directly from factors
% means=reshape(u*fmeans',nx,ny,nclasses);

% Compute full-sized means
means=single(zeros(nx0,ny0,nclasses));
for i=1:nclasses
    means(:,:,i)=squeeze(mean(img(:,:,indices==i),3));
end;

nmembers=hist(indices,1:nclasses);
fmeans=fmeans';  % transpose for returned value.

% figure(1);
% ImagicDisplay(means);
%
% figure(2); clf;
% bar(nmembers);
