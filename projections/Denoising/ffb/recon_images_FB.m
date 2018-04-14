% RECON_IMAGES_FB Reconstruct images from Fourier-Bessel coefficients
%
% Usage
%    recon = recon_images_FB(c, R, L0, Coeff, n_min, n_max);
%
% Input
%    c: The band limit of the basis in frequency.
%    R: The radius of the disk on which the basis is concentrated.
%    L0: The size of the images.
%    Coeff: The Fourier-Bessel coefficients in the same format as FBcoeff_nfft.
%    n_min, n_max: The range of images indices to reconstruct.
%
% Output
%    recon: An array of size L0-by-L0-by-(n_max-n_min+1) containing the
%       reconstructed images.

% Written by Zhizhen Zhao - 10/15
% Reformatted and documented by Joakim Anden - 2018-Apr-13

function recon = recon_images_FB(c, R, L0, Coeff, n_min, n_max)
    L = 2*R;

    fn = IFT_FB(R, c);
    max_ang_freqs = size(fn, 1)-1;
    origin = floor(L0/2)+1;

    tmp1 = fn{1};
    tmp1 = reshape(tmp1, L^2, size(tmp1, 3));

    recon = zeros(L0, L0, n_max-n_min+1);
    for K = n_min:n_max
        tmp2 = tmp1*Coeff{1}(:,K);
        tmp2 = reshape(tmp2, L, L);
        I = real(tmp2);
        for k = 1:max_ang_freqs
            tmp = fn{k+1};
            tmp = reshape(tmp, L^2, size(tmp, 3));
            tmp2_pos = tmp*Coeff{k+1}(:,K);
            tmp_2 = 2*real(tmp2_pos);
            I = I + reshape(tmp_2, L, L);
        end

        test = zeros(L0);
        test(origin-R:origin+R-1,origin-R:origin+R-1) = real(I);
        recon(:,:,K-n_min+1) = test;
    end
end
