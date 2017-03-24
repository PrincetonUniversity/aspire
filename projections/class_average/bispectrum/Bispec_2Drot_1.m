function [ Coeff_b,  toc_bispec ] = Bispec_2Drot( Coeff, Freqs )
%Description:
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
%       toc_bispec:
%           Timing for bispectrum computation.
%
%   Zhizhen Zhao Aug 2013
%   Tejal Updated April 2016
tic_bispec=tic;
alpha=1/3; %modify the amplitude for each component.
n_im = size(Coeff, 2);
Coeff_norm=(abs(Coeff(Freqs~=0, :))).^(alpha);
Coeff_norm=log(Coeff_norm + 1e-15);
%Adding epsilon to zero coefficients before taking log, Tejal April 18, 16
check=isinf(Coeff_norm);
assert(max(check(:))~=1);
clear check

Phase=Coeff(Freqs~=0, :)./abs(Coeff(Freqs~=0, :));
check=isnan(Phase);
if(max(check(:))==1)
	Phase=atan2(imag(Phase), real(Phase));
	Phase(isnan(Phase))=0; % Will happen if some coeffs are exactly zero
else
	assert(max(check(:))~=1);
	clear check
	Phase=atan2(imag(Phase), real(Phase));
end

clear Coeff
[O1, O2] = bispec_Operator_1(Freqs(Freqs~=0));
fprintf('\nLength of the bispectrum is %d \n', size(O1, 1));

M=exp(O1*Coeff_norm+sqrt(-1)*O2*Phase);
ncomp = min(min(size(M,1), size(M,2)),200);

check=isnan(M);
assert(max(check(:))~=1);
clear check

%dimensionality reduction using PCA
[ U, S, V ]=pca_Y(M, ncomp);
%bar(diag(S));
fprintf('\nFinished PCA');
Coeff_b = S*V';
%Normalize 

for i=1:size(Coeff_b, 2);
     Coeff_b(:, i)=Coeff_b(:, i)/norm(Coeff_b(:, i));
end;
toc_bispec=toc(tic_bispec);
end

