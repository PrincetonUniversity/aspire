function [ coeff_pos_k ]= FBcoeff_nfft(data, R, basis, sample_points, num_pool)
%Description:
%This code use nufft to compute Fourier Bessel expansion coefficients for
%k>=0.
%Input: 
%   data: images data LxLxn
%   R: compact support size in real space
%   basis: Bessel basis computed in precomp.m
%   sample_points: GQ points and weights, computed in precomp.
%Output:
%   coeff_pos_k: Fourier bessel expansion coefficients in cell structure.
%   number of cells = k_max + 1. coeff_pos_k{i} contain coefficients for k
%   = i-1.
%
%Zhizhen Zhao 09/2015

%num_pool = maxNumCompThreads;
%num_pool = 30;
%n_im = size(data, 3);
L = size(data{1}, 1); %initial image size
orig = floor(L/2)+1;
nL = 2*R;  % final image after cropping the image
%Putting images into batches

%Read input data
Phi_ns = basis.Phi_ns;
ang_freqs = basis.ang_freqs;
n_theta = basis.n_theta;
clear basis;

max_ang_freqs = max(ang_freqs);

w = sample_points.w;
r = sample_points.r;
w = r.*w;

[ freqs ] = pft_freqs(r, n_theta);
Precomp.n_theta = n_theta;
Precomp.n_r = length(r);
Precomp.nL = nL;
Precomp.freqs = freqs;
scale = 2*pi/(n_theta);

%Evaluate expansion coefficients

coeff_pos_k = cell(max_ang_freqs+1, 1);
pos_k = cell(num_pool, 1);

%poolreopen(num_pool);
%disp('Setting OMP_NUM_THREADS to 1 for nufft')
%setenv('OMP_NUM_THREADS', '1');

parfor i = 1:num_pool
    tmp = data{i};
    tmp = tmp(orig-R:orig+R-1, orig-R : orig+R-1, :);
    tmp2 = cryo_pft_nfft(tmp, Precomp);
    pf_f = scale*fft(tmp2, [], 2); %1D FFT on concentric rings
    pos_k{i} = pf_f(:, 1:max(ang_freqs)+1, :);
end;

pos_k = cat(3, pos_k{:});

clear pf_f tmp tmp2;

for i = 1:max_ang_freqs+1
    tmp = squeeze(pos_k(:, i, :));
    coeff_pos_k{i} = Phi_ns{i}.'*diag(w)*tmp;
end;

end
