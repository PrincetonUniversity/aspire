function [ Coeff_b, Coeff_b_r, toc_bispec ] = Bispec_2Drot( Coeff, Freqs )
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
Coeff_norm=(abs(Coeff(Freqs~=0, :))).^(alpha);
Coeff_norm=log(Coeff_norm);
check=isinf(Coeff_norm);
assert(max(check(:))~=1);
clear check
Phase=Coeff(Freqs~=0, :)./abs(Coeff(Freqs~=0, :));
Phase=atan2(imag(Phase), real(Phase));
[O1, O2]=bispec_Operator(Freqs(Freqs~=0));
M=exp(O1*Coeff_norm+sqrt(-1)*O2*Phase);
clear Coeff_norm Phase O1 O2

%% PCA the bispectrum

[U, S, V]=pca_Y(M, min([200 size(M)]));
Coeff_b=S*V';       %reduced rotationally invariant features
Coeff_b_r=U'*conj(M);      %reduced rotationally invariant features for reflected images.
clear M U S V

%figure; bar(diag(S));
for i=1:size(Coeff_b, 2);
    Coeff_b(:, i)=Coeff_b(:, i)/norm(Coeff_b(:, i));
    Coeff_b_r(:, i)=Coeff_b_r(:, i)/norm(Coeff_b_r(:, i));
end;
toc_bispec=toc(tic_bispec);
end

