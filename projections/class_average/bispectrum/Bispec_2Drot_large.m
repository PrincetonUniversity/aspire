function [ Coeff_b, Coeff_b_r, toc_bispec ] = Bispec_2Drot_large( Coeff, Freqs )
%This function computes rotationally invariant bispectral-like features for 2D images.
%   Input: 
%       Coeff:
%           Truncated expansion coefficients of 2D images on Fourier Bessel
%           steerable PCA basis.
%       Freqs:
%           The angular frequencies associated with each component.
% 
%   Output:
%       Coeff_b:
%           The invariant feature for original images.
%       Coeff_b_r:
%           The invariant feature for the reflected images.
%       toc_bispec:
%           Timing for bispectrum computation.
%
%   Zhizhen Zhao Aug 2013

tic_bispec=tic;
alpha=1/3; %modify the amplitude for each component.
N = size(Coeff, 2);
Coeff_norm=(abs(Coeff(Freqs~=0, :))).^(alpha);
Coeff_norm=log(Coeff_norm);
check=isinf(Coeff_norm);
assert(max(check(:))~=1);
clear check
Phase=Coeff(Freqs~=0, :)./abs(Coeff(Freqs~=0, :));
Phase=atan2(imag(Phase), real(Phase));
[O1, O2] = bispec_Operator_1(Freqs(Freqs~=0));

%For large data, use the top 10000 to approximate the low dimensional
%projecion

BatchSize = 10000;
nB = ceil(N/BatchSize);
M=exp(O1*Coeff_norm(:, 1:min(N, BatchSize))+sqrt(-1)*O2*Phase(:, 1:min(N, BatchSize)));

%% PCA the bispectrum
[ U, ~, ~ ] = pca_Y( M, min([200 size(M)]) );
Coeff_b = zeros( min([200 size(M)]), N);
Coeff_b_r = zeros( min([200 size(M)]), N);
for i = 1:nB
    M = exp(O1*Coeff_norm(:, (i-1)*BatchSize+1: min(i*BatchSize, N)) + sqrt(-1)*O2*Phase(:, (i-1)*BatchSize+1: min(i*BatchSize, N)));
    Coeff_b(:, (i-1)*BatchSize+1: min(i*BatchSize, N)) = U'*M;       %reduced rotationally invariant features
    Coeff_b_r(:, (i-1)*BatchSize+1: min(i*BatchSize, N)) = U'*conj(M);      %reduced rotationally invariant features for reflected images.
end;

for i=1:size(Coeff_b, 2)
    Coeff_b(:, i) = Coeff_b(:, i) / norm(Coeff_b(:, i));
    Coeff_b_r(:, i) = Coeff_b_r(:, i) / norm(Coeff_b_r(:, i));
end;
toc_bispec=toc(tic_bispec);

end

