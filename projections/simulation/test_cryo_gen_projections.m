% Test the function cryo_add shifts.
%
% Yoel Shkolnisky, September 2013.

K=10;
SNR=1;
max_shift=5;
shift_step=1;
n=65;

tic;
[clean_projections_2, noisy_projections_2, shifts, q] = ...
    cryo_gen_projections(n,K,SNR);
toc

