% CRYO_PERIODOGRAM_COVARIANCE Calculate covariance of periodograms
%
% Usage
%    x_per_cov = cryo_periodogram_covariance(x_per);
%
% Input
%    x_per: An array of periodograms of size L1-by-L2-by-n obtained from
%       cryo_periodogram.
%
% Output
%    x_per_cov: An image matrix of size L1-by-L2-by-L1-by-L2 consisting of
%       the covariances between the various frequencies in the periodograms.

function x_per_cov = cryo_periodogram_covariance(x_per)
    n = size(x_per, 3);

    L = [size(x_per, 1) size(x_per, 2)];

    x_per_half = fourier_full_to_half(x_per, L);

    C = x_per_half*x_per_half'/n;

    x_per_mean = cryo_periodogram_mean(x_per);
    x_per_mean = fourier_full_to_half(x_per_mean, L);

    [mask, real_mask] = positive_fourier_mask(L);

    mult = ones(sum(mask(:))) + diag(1+real_mask(mask));

    x_per_cov = C./mult-x_per_mean*x_per_mean';

    x_per_cov = fourier_half_to_full(x_per_cov, L);
    x_per_cov = permute(x_per_cov, [3 1:2]);
    x_per_cov = fourier_half_to_full(x_per_cov, L);
    x_per_cov = permute(x_per_cov, [3:4 1:2]);
end
