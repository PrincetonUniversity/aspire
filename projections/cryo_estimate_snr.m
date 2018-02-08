function [snr,var_s,var_n]=cryo_estimate_snr(projs,prewhiten)
% CRYO_ESTIMATE_SNR     Estimate SNR of projections stack
%
% [snr,var_s,var_n]=cryo_estimate_snr(projs)
%   Estimate the snr of a stack of projections by using corner pixels. If
%   each projections in the stack is of size pxp, then only pixels whose
%   radius is larger than p/2 are used for estimating the parameters of the
%   noise. snr is defined as var_s/var_n, where var_s is the variance of
%   the signal and var_n is the variance of the noise.
%   Noise in the proejctions is assumed to be white.
%
% [snr,var_s,var_n]=cryo_estimate_snr(projs,prewhiten)
%   If prewhiten is non-zero, then prewhiten the projections before
%   estimating the parameters of the noise.
%   
%
% IMPORTANT: All images must be normalized to have the same noise mean and
% variance, using, e.g., cryo_normalize_background
%
% Yoel Shkolnisky, August 2015.


% prewhiten if needed
if nargin==2
    if prewhiten
        psd = cryo_noise_estimation(projs);
        prewhitened_projs = cryo_prewhiten(projs, psd);
        projs=prewhitened_projs;
        clear prewhitened_projs;
    end
end

% Masking radius
p = size(projs, 1); % Image size)
radius_of_mask = floor(p/2)-1; % 

% Compute indices of noise samples
I = cart2rad(p);
mask = zeros(p);
mask(I<radius_of_mask) = -1;  % Create a mask.
noise_idx = find(mask~=-1);   % Indices of noise pixels

% Compute noise variance from all projections
%sum_noise=0;    % Sum of all noise pixels;
%sum_noise_sq=0; % Sum of noise pixels squared
total_noise_samples=size(projs,3)*numel(noise_idx); % Number of noise sample

% Calculate the mean of the noise and the signal+noise.
noise_means=zeros(size(projs,3),1,class(projs));
projs_means=zeros(size(projs,3),1,class(projs));
parfor k=1:size(projs,3)
    current_proj=projs(:,:,k);
    noise_vals=current_proj(noise_idx);
%    sum_noise=sum_noise+sum(noise_vals);
    noise_means(k)=mean(noise_vals);
    projs_means(k)=mean(current_proj(:));
end
%mean_noise=sum_noise/total_noise_samples; % Noise mean
%noise_mean=mean(noise_means);

% Calculate the sample variance of the noise and the projections.
% Due to the normalization of the projections, there is no need to estimate
% the variance of the noise, as it is equal to 1. However, I estimate it to
% verify that nothing went wrong.
noise_sumsq=zeros(size(projs,3),1,class(projs));
projs_sumsq=zeros(size(projs,3),1,class(projs));
parfor k=1:size(projs,3)
    current_proj=projs(:,:,k);
    noise_vals=current_proj(noise_idx);
%    sum_noise_sq=sum_noise_sq+sum((abs(noise_vals-mean_noise)).^2);
    noise_sumsq(k)=sum((abs(noise_vals-noise_means(k))).^2);
    
    % Note that it is incorrect to subtract the mean image, since then the
    % code won't work for a stack consisting of multiple copies of the same
    % image plus noise.
    projs_sumsq(k)=sum((abs(current_proj(:)-projs_means(k))).^2);
end
var_n=sum(noise_sumsq)/(total_noise_samples-1);
var_splusn=sum(projs_sumsq)/(numel(projs)-1);
%var_n=sum_noise_sq/(total_noise_samples-1); % Noise variance

%mean_proj = mean(projs, 3);
%mean_proj=cryo_radial_average2d(mean_proj);
%projs = bsxfun(@minus, projs, mean_proj);

%var_s=mean(projs(:).^2)-var_n; % Signal variance
var_s=var_splusn-var_n;
snr=var_s/var_n;

