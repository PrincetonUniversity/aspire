function [ coeff_pos_k ]= FBcoeff_sphere(data_image, basis, sample_points)
%Description:
%This code use nufft to compute Fourier Bessel expansion coefficients for
%k>=0.
%Input: 
%   data: images data n_r*L^2
%   basis: Bessel basis computed in precomp.m
%   sample_points: GQ points and weights, computed in precomp.
%Output:
%   coeff_pos_k: Fourier bessel expansion coefficients in cell structure.
%   number of cells = k_max + 1. coeff_pos_k{i} contain coefficients for k
%   = i-1. In coeff_pos_k{i,j}, j stands for the lm number in the spherical
%   harmonics expansion.
%
%Zhizhen Zhao 09/2015
% Yuan Liu 11/2016

Phi_ns = basis.Phi_ns;
w = sample_points.w;
r = sample_points.r;
w = r.^2.*w;

%Evaluate expansion coefficients
Lmax = length(Phi_ns);
coeff_pos_k = cell(Lmax, 1);

for l = 1:Lmax
    m = (l-1)^2+(1:2*l-1);
    coeff_pos_k{l,1} = Phi_ns{l}.'*diag(w)*data_image(:,m); 
end

end
