function [ Coeff_b, Coeff_b_r, toc_bispec ] = Bispec_2Drot_large( Coeff, Freqs, eigval )
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
eigval = eigval(Freqs~=0);
[O1, O2] = bispec_Operator_1(Freqs(Freqs~=0));
clear Coeff;
%disp('cleared coefficients');
%generate the variance of the bispectral coefficients
N = 4000;
M = exp(O1*log(eigval.^(alpha)));
pM = M/sum(M);
x = rand(length(M), 1);
M_id = find(x<N*pM); 
log_message('Number of bispectrum components is %d',length(M_id));
%[ ~, M_id ] = sort(M, 'descend');
O1 = O1(M_id, :);
O2 = O2(M_id, :);

M=exp(O1*Coeff_norm+sqrt(-1)*O2*Phase);

%% svd of the reduced bispectrum
[ U, S, V ] = pca_Y( M, 300 );
%disp('PCA done');
Coeff_b = S*V';
Coeff_b_r = U'*conj(M);

for i=1:size(Coeff_b, 2)
    Coeff_b(:, i) = Coeff_b(:, i) / norm(Coeff_b(:, i));
    Coeff_b_r(:, i) = Coeff_b_r(:, i) / norm(Coeff_b_r(:, i));
end;
toc_bispec=toc(tic_bispec);

end

