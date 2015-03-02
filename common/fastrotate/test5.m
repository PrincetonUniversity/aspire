function test5
% Test the function fastrotate with precomputation.
%
% Rotate a Gaussian 360 degrees and check we get back the same Gaussian.
% The error should be around 10^-13.
% This file cannot be a script since we use an inner function.
%
% Yoel Shkolnisky, January 2011.

% Even test
N=128;
runtest(N);
disp('Press any key to continue');
pause;

% Even test
N=129;
runtest(N);

% Run a test with an NxN image
function runtest(N)
c=N/2;
[X,Y]=meshgrid(0:N-1,0:N-1);
X=X-c;
Y=Y-c;
sigmax=N/20;
sigmay=N/30;
im=exp(-(X.^2/(2*sigmax.^2)+Y.^2/(2*sigmay.^2)));

subplot(1,3,1);
imagesc(im);
axis image;
axis off;

M=fastrotateprecomp(N,N,1);
rim=im;
for k=1:360;
    rim=fastrotate(rim,[],M);
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
