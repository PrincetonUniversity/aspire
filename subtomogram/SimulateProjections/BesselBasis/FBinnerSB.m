function [ innerP ]= FBinnerSB(basis_fb, basis, sample_points)
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


%Read input data
Phi_fb = basis_fb.Phi_ns;
Phi_ns = basis.Phi_ns;
clear basis;

w = sample_points.w;
r = sample_points.r;
w = r.*w;

%Evaluate inner product
L = size(Phi_ns,1);
innerP = cell(L,L);

for i = 1:size(Phi_ns,1)
    temp = cell(1,L);
    for j = 1:i
        temp{1,j} = Phi_ns{i}.'*diag(w)*Phi_fb{j};
    end
    innerP(i,:) = temp;
end;


end
