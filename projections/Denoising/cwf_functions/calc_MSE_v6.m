function [mse, mse_all] = calc_MSE_v6(im, im_true, r_max)

% Calculate MSE (squared Frobenius norm) inside a masked region of radius r_max = floor(L/2)
% im: Estimated image stack
% im_true: True image stack
% mse is the average MSE across all images
% Tejal Bhamre: April 2015

L=size(im,1);
nim=size(im,3);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
% r_max=floor(L/2);

for i=1:nim
    temp=zeros(L,L);
    I=im_true(:,:,i);
    temp(r<=r_max)=I(r<=r_max);
    I_masked(:,:,i)=temp;
    den_im=im(:,:,i);
    temp(r<=r_max)=den_im(r<=r_max);
    den_masked(:,:,i)=temp;
end

real_ims=real((reshape(I_masked,L*L,nim)));
denoised_ims=real(reshape(den_masked,L*L,nim));
err=(real_ims-denoised_ims);
%mse=mean(sum(err.^2,1));

ims_sqnorm=sum(real_ims.^2,1);
mse=mean(sum(err.^2,1)./ims_sqnorm); % Relative MSE, MSE/norm of true image squared
mse_all=sum(err.^2,1)./ims_sqnorm;
