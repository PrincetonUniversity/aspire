function [ mean_image, denoised] = denoise_images_analytical_fix(U, fn, mean_coeff, Coeff, L0, R, n_max)
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

max_ang_freqs = size(U, 1)-1; %Should be the maximum angular frequency
L = 2*R;
%Computes eigen images, need output from IFT_FB.m.
eig_im = cell(max_ang_freqs+1, 1);
for k = 1:max_ang_freqs+1
    if size(U{k},2)~=0 
        tmp = fn{k};
	tmp = reshape(tmp, L^2, size(tmp, 3));
        tmp2 = tmp*U{k};
        eig_im{k} = reshape(tmp2, L, L, size(tmp2, 2));
    end;
end;

%Original image size
origin = floor(L0/2) + 1;

tmp = fn{1};
tmp = reshape(tmp, L^2, size(tmp, 3));
mean_Im = reshape(tmp*mean_coeff, L, L);
mean_Im = real(mean_Im);

if(size(eig_im{1},2)~=0)
tmp1 = eig_im{1};
tmp1 = reshape(tmp1, L^2, size(tmp1, 3));
end


denoised = zeros(L0, L0, n_max);
for K = 1:n_max

%% Added this condition to fix the case when k=0 component is not picked
%Tejal

if(size(eig_im{1},2)~=0)  
tmp2 = tmp1*Coeff{1}(:, K);
tmp2 = reshape(tmp2, L, L);
I = mean_Im + real(tmp2);
else
I = mean_Im;
end
    
for k = 1:max_ang_freqs
    if size(U{k+1},2)~=0 
        tmp = eig_im{k+1};
        tmp = reshape(tmp, L^2, size(tmp, 3));
        tmp2_pos = tmp*Coeff{k+1}(:, K);
	tmp_2 = 2*real(tmp2_pos);
        I = I + reshape(tmp_2, L, L);
    end;
end;

test = zeros(L0);
test(origin-R:origin+R-1, origin-R:origin+R-1) = real(I);
denoised(:, :, K) = test;

end;

mean_image = zeros(L0);
mean_image(origin-R:origin+R-1, origin-R:origin+R-1) = mean_Im;

