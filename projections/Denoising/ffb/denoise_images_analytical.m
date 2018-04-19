% DENOISE_IMAGES_ANALYTICAL Mean and denoised images from sPCA coefficients
%
% Usage
%    [mean_image, denoised] = denoise_images_analytical(U, fn, mean_coeff, ...
%        coeff, L, R, N);
%
% Input
%    U: A cell array containing the eigenvectors of the covariance blocks
%       C^{(1)} through C^{(k_max)}. The kth element contains the eigenvectors
%       C^{(k)}, arranged in columns.
%    fn: The inverse Fourier transform of the Fourier-Bessel basis vectors. This
%       is obtained using the function IFT_FB.
%    mean_coeff: The mean Fourier-Bessel coefficients for k = 0.
%    coeff: The sPCA coefficients obtained from Wiener filtering.
%    L0: The desired size of the output images.
%    N: The number of coefficients to process.
%
% Output
%    mean_image: The mean image in the spatial domain in a grid of size
%       L0-by-L0.
%    denoised: The denoised images, arranged in an array of size L0-by-L0-by-N.

% Written by Zhizhen Zhao - 10/15
% Updated by Tejal - 08/17
% Reformatted and documented by Joakim Andén - 2018-Apr-19

function [mean_image, denoised] = denoise_images_analytical(U, fn, mean_coeff, coeff, L, R, N)
    max_ang_freqs = size(U, 1)-1;

    % Computes eigenimages, needs output from IFT_FB.m.
    eig_im = cell(max_ang_freqs+1, 1);
    for k = 1:max_ang_freqs+1
        if size(U{k}, 2) ~= 0
            tmp = fn{k};
            tmp = reshape(tmp, (2*R)^2, size(tmp, 3));
            tmp2 = tmp*U{k};
            eig_im{k} = reshape(tmp2, 2*R, 2*R, size(tmp2, 2));
        end
    end

    origin = floor(L/2) + 1;
    tmp1 = fn{1};
    tmp1 = reshape(tmp1, (2*R)^2, size(tmp1, 3));
    mean_Im = reshape(real(tmp1*mean_coeff), 2*R, 2*R);

    tmp1 = eig_im{1};
    tmp1 = reshape(tmp1, (2*R)^2, size(tmp1, 3));
    tmp2 = tmp1*coeff{1}(:,1:N);
    tmp2 = reshape(tmp2, 2*R, 2*R, N);
    I = mean_Im + real(tmp2);

    for k = 1:max_ang_freqs
        if size(coeff{k+1},1) ~= 0
            tmp = eig_im{k+1};
            tmp = reshape(tmp, (2*R)^2, size(tmp, 3));
            tmp2_pos = tmp*coeff{k+1}(:,1:N);
            tmp_2 = 2*real(tmp2_pos);
            I = I + reshape(tmp_2, 2*R, 2*R, N);
        end
    end
    denoised = zeros(L, L, N);
    denoised(origin-R:origin+R-1,origin-R:origin+R-1,:) = real(I);

    % Estimated mean images
    mean_image = zeros(L);
    mean_image(origin-R:origin+R-1,origin-R:origin+R-1) = mean_Im;
end
