function [ timing, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(data, R, noise_variance, basis,  sample_points, num_pool)
%Description:
% Computes Fourier-Bessel expansion coefficients 'coeff', and filtered steerable PCA expansion coefficients 'sPCA_coeff'.
%Input:
%	data: image dataset.
%	c: band limit
%	R: compact support radius
%	noise_variance: estimated noise variance
%Output
%	timing: toc_FFBsPCA: the whole time spent on computing fast steerable PCA. toc_sPCA: time for steerable PCA. toc_FBcoeff: time for computing Fourier-Bessel expansion coefficients
%	coeff: Fourier-Bessel expansion coefficients
%	mean_coeff: mean of the Fourier-Bessel expansion coefficients
%	sPCA_coeff: steerable PCA expansion coefficients
%	U: eigenvectors of C^{(k)}'s
%	D: eigenvalues of C^{(k)}'s. 
%Update: 10/15 Zhizhen Zhao

%n_r = ceil(4*c*R);
n = size(data, 3);
data2 = cell(num_pool, 1);
nb = floor(n/num_pool);
remain = n - nb*num_pool;
for i = 1:remain
    data2{i} = data(:, :,(nb+1)*(i-1)+1: (nb+1)*i);
end;
count = (nb+1)*remain;
for i = remain+1:num_pool
    data2{i} = data(:, :, count + (i-remain-1)*nb+1: count + (i-remain)*nb);
end; 
clear data;
tic_start = tic;
[coeff ]= FBcoeff_nfft(data2, R, basis, sample_points, num_pool);
disp('Finished computing Fourier Bessel coefficients')
toc_FBcoeff = toc(tic_start);
tic_start2 = tic;
disp('Steerable PCA: computing compressed coefficients')
[ U, D, sPCA_coeff, mean_coeff ] = sPCA_whole( coeff, basis, noise_variance);
toc_sPCA = toc(tic_start2);
toc_FFBsPCA = toc_sPCA + toc_FBcoeff;

timing.toc_FBcoeff = toc_FBcoeff;
timing.toc_sPCA = toc_sPCA;
timing.toc_FFBsPCA = toc_FFBsPCA;
end
