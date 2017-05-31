function [noisy_real, noise_v]=addnoise_v6(signal, SNR)

% Add gaussian noise to image stack
% INPUT:
% signal: Images with CTF and/or envelope in fourier domain
% SNR: Desired SNR of noisy images
% OUTPUT:
% noise_v: Noise variance
% noisy_real: Stack of noisy images
% Tejal Bhamre: updated Oct 22/15

rng('default')
noise_std=sqrt(var(signal(:))/SNR);   % noise_std is standard dev
noisy_real=signal+randn(size(signal))*noise_std;
noise_v=noise_std.^2;
