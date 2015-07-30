% Test the function register_translations_2d
%
% Yoel Shkolnisky, January 2014.

%% Load an image
data=load('projs');
im1=data.projs(:,:,1);

% Make the image even-sized.
% Note that for even-sized image there are small errors when the shifts are
% not integral pixels (of the order of 10^-5).
%im1=im1(1:end-1,1:end-1);

%% Create a shifted image
% Integral shift, without calling reshift_image
%im2=im1([4:end 1:3],[3:end 1:2]);
%dx=[3 2];

dx=[1.543,3.777];  % This is the shift we need to recover.
%dx=[2,0];
im2=reshift_image(im1,dx);  % First shift paramter is the x coordinate. 

%% Add noise to the images.
sigma=0; % If sigma=0 (no noise), then the shift dx should be recovered very accurately.
im1=im1./norm(im1(:));
im2=im2./norm(im2(:));

rng('default')
rng(111);
n1=sigma*randn(size(im1))./sqrt(numel(im1(:)));
n2=sigma*randn(size(im2))./sqrt(numel(im2(:)));
im1=im1+n1;
im2=im2+n2;

% Show the two images to resiter
h1=figure;
subplot(1,3,1);
imagesc(im1); axis image; axis off; colormap(gray);
title('im1');
subplot(1,3,2);
imagesc(im2); axis image; axis off; colormap(gray);
title(im2);

%% Register 
estdx=register_translations_2d(im1,im2,dx);
figure(h1);
subplot(1,3,3);
imagesc(reshift_image(im2,estdx)-im1); axis image; axis off;  colorbar;
%title('Difference between im1 and im2 after registration');

% % % % Register with filtering.
% % % % No point in filtering. Phase correlation eliminates any linear filter.
% % % im1=GaussFilt(im1,0.1);
% % % im2=GaussFilt(im2,0.1);
% % % h2=figure;
% % % subplot(1,3,1)
% % % imagesc(im1); axis image; colormap(gray);
% % % subplot(1,3,2)
% % % imagesc(im2); axis image; colormap(gray);
% % % 
% % % estdx=register_translations_2d(im1,im2,dx);
% % % 
% % % figure(h2);
% % % subplot(1,3,3)
% % % imagesc(reshift_image(im2,estdx)-im1);axis image; colorbar;