function [mean_image, denoised] = denoise_images_analytical(U, fn, mean_coeff, coeff, L, R, N)
% Description
% This function gives mean and denoised images in real space
%
% Input: 
%	U: eigenvectors of C^{(k)}, saved in cell array
% 	fn: inverse Fourier transform of the Fourier- Bessel function in real space on a Cartesian grid.
%	mean_coeff: mean Fourier-Bessel expansion coefficients for k=0
%	Coeff: sPCA expansion coefficients (Wienter filtered).
%	L0: original image size
%	n_max: number of images you want to denoise
%	 
% Output:
%	mean_image: mean image in real space on a Cartesian grid of size L0 x L0
%	denoised: denoised images in real space on a Cartesian grid L0 x L0
% Update 10/15 Zhizhen Zhao
% Updated Tejal
% Updated 08/17
% See Fast sPCA paper

max_ang_freqs = size(U, 1)-1; %Should be the maximum angular frequency
%Computes eigen images, need output from IFT_FB.m.
eig_im = cell(max_ang_freqs+1, 1);
for k = 1:max_ang_freqs+1
    if size(U{k},2)~=0 
        tmp = fn{k};
	    tmp = reshape(tmp, (2*R)^2, size(tmp, 3));
        tmp2 = tmp*U{k};
        eig_im{k} = reshape(tmp2, 2*R, 2*R, size(tmp2, 2));
    end;
end;

origin = floor(L/2) + 1;
tmp1 = fn{1};
tmp1 = reshape(tmp1, (2*R)^2, size(tmp1, 3));
mean_Im = reshape(real(tmp1*mean_coeff), 2*R, 2*R);

tmp1 = eig_im{1};
tmp1 = reshape(tmp1, (2*R)^2, size(tmp1, 3));
tmp2 = tmp1*coeff{1}(:,  1:N );
tmp2 = reshape(tmp2, 2*R, 2*R, N);
I = mean_Im + real(tmp2);

for k = 1:max_ang_freqs
    if size(coeff{k+1},1)~=0
        tmp = eig_im{k+1};
        tmp = reshape(tmp, (2*R)^2, size(tmp, 3));
        tmp2_pos = tmp*coeff{k+1}(:, 1 : N);
        tmp_2 = 2*real(tmp2_pos);
        I = I + reshape(tmp_2, 2*R, 2*R, N);
    end;
end;
denoised = zeros( L, L, N );
denoised(origin-R:origin+R-1, origin-R:origin+R-1, :) = real(I);

%estimated mean images
mean_image = zeros(L);
mean_image(origin-R:origin+R-1, origin-R:origin+R-1) = mean_Im;
