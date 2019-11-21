function [recon] = recon_images_FBm_2(R, L0, Coeff, n_min, n_max, fn2, r0)
% Description
% This function reconstruct image in real space from Fourier-Bessel expansion coefficients
% Input:
%	c: band limit
%	R: compact support radius R
%	L0: original size of the images
%	Coeff: Fourier-Bessel expansion coefficients in cell structure
  %     n_min: Starting image
%	n_max: End image
% Output:
%	recon: reconstructed images in real domain
% Update: 10/15 Zhizhen Zhao 

%provide bandlimit c and compact support size R
%provide expansion coefficients 'Coeff'
%L = 2*R+1;

%max_ang_freqs = size(Coeff, 1)-1; %find k_max Yuan:changed from fn to Coeff
n_im = n_max-n_min+1;

recon2 = fn2*Coeff;
recon = zeros(L0*L0,n_im);
recon(r0(:)<=R,:) = real(recon2);
recon = reshape(recon, L0, L0, n_im);
recon = recon*(2*pi)^(1/2);

end

