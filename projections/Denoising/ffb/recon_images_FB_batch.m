function [recon] = recon_images_FB_batch(c, R, L0, batchid)
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
L = 2*R;
%Computes eigen images, need output from IFT_FB.m.
[ fn ] = IFT_FB(R, c);
max_ang_freqs = size(fn, 1)-1; %find k_max
origin = floor(L0/2) + 1;

tmp1 = fn{1};
tmp1 = reshape(tmp1, L^2, size(tmp1, 3));

filename=fullfile('/scratch/tbhamre/cwf_batch/', sprintf('batch_den_coeff_pos%d',batchid));
load(filename);
n_min=1;
n_max=size(den_coeff_new{1},2);

recon = zeros(L0, L0, n_max-n_min+1);
for K = n_min:n_max
    tmp2 = tmp1*den_coeff_new{1}(:, K);
    tmp2 = reshape(tmp2, L, L);
    I = real(tmp2);
    for k = 1:max_ang_freqs
        tmp = fn{k+1};
        tmp = reshape(tmp, L^2, size(tmp, 3));
        tmp2_pos = tmp*den_coeff_new{k+1}(:, K);
        tmp_2 = 2*real(tmp2_pos);
        I = I + reshape(tmp_2, L, L);
    end;

    test = zeros(L0);
    test( origin-R : origin+R - 1, origin-R : origin+R-1) = real(I);
    recon(:, :, K-n_min+1) = test;
%    sprintf('Done reconstructing %dth image',K)

end;

