% Test Yaroslavsy's rotation code.
%
% Same as test1, but with a function that is NOT smooth.
% This demosntrates the artifacts for images with sharp edges.
%
% Yoel Shkolnisky, January 2011.

N=128;
c=N/2;
im=zeros(N);
%im(c-N/4+1:c+N/4,c-N/4+1:c+N/4)=1;
im(c-N/8+1:c+N/8,c-N/8+1:c+N/8)=1;

subplot(1,3,1);
imagesc(im);
axis image;
axis off;

rim=im;
for k=1:360;
    rim=fastrotate_yaro(rim,1);
    subplot(1,3,2);
    imagesc(rim);
    axis image;
    axis off;
    pause(0.01);
end

subplot(1,3,3);
imagesc(rim-im);
axis image;
axis off;
colorbar;
fprintf('err=%e\n',norm(rim(:)-im(:))/norm(im(:))); 