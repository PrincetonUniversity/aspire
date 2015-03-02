% Test Yaroslavsy's rotation code.
%
% Rotate and narrow Gaussian and check we get back the same Gaussian.
% Note that if the Gaussian is narrow, the error is very small.
%
% Yoel Shkolnisky, January 2011.

N=128;
c=N/2;
[X,Y]=meshgrid(0:N-1,0:N-1);
X=X-c;
Y=Y-c;
sigmax=N/10;
sigmay=N/20;
im=exp(-(X.^2/(2*sigmax.^2)+Y.^2/(2*sigmay.^2)));

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
disp(norm(rim(:)-im(:))/norm(im(:))); 
    
    