function [ fn ] = IFT_FB(R, c)

%%% Description
%The function computes the inverse Fourier transform of Fourier-Bessel Basis
%   This evaluate Bessel radial functions for an image of size
%   Table bessel.mat gives the R_ns (3rd column) for n (1st column) and s (2nd column)
% Input:  
%       1. c: band limit
%       2. R: compact support radius in real domain
% Output: 
%       1. fn: IFT of FB basis, stored in cell structure. cell number corresponds to the angular index.     
% Zhizhen Zhao 04/2015

[ x, y ] = meshgrid(-R:R-1, -R:R-1);
r = sqrt(x.^2 + y.^2);
r0 = r;
r = r(:);
theta = atan2(y, x);
theta = theta(:);
theta = theta(r(:)<=R);
r = r(r(:)<=R);

load bessel.mat
bessel = bessel(bessel(:, 4)<= 2*pi*c*R, :);
k_max = max(bessel(:, 1));
fn = cell(k_max+1,1);
for i = 1:k_max+1
    bessel_k = bessel(bessel(:, 1) == i-1, :);
    l = size(bessel_k, 1);
    tmp2 = zeros(2*R);
    tmp3 = zeros(2*R, 2*R, l);
    for lr = 1:l
    	Jk = besselj(i-1, 2*pi*c*r(:));
    	Rkq = bessel_k(lr, 3);
    	f_r = 2*c*sqrt(pi)*(-sqrt(-1))^(i-1)*(-1)^lr*Rkq*Jk./((2*pi*c*r).^2 - Rkq^2);
    	f_theta = exp(sqrt(-1)*(i-1)*theta(:));
    	tmp = f_r.*f_theta;
        tmp2(r0<=R) = tmp;
    	tmp3(:, :, lr) = tmp2;
    end;
    fn{i} = tmp3;
end;

%for i = 1:k_max
%    bessel_k = bessel(bessel(:, 1)==i, :);
%    l = size(bessel_k, 1);
%    tmp2 = zeros(2*R);
%    tmp3 = zeros(2*R, 2*R, l);
%    for lr = 1:l
%    	Jk = besselj(i, 2*pi*c*r(:));
%    	Rkq = bessel_k(lr, 3);
%    	f_r = 2*c*sqrt(pi)*(sqrt(-1))^(i)*(-1)^lr*Rkq*Jk./((2*pi*c*r).^2 - Rkq^2 );
%        f_theta = exp(-sqrt(-1)*(i)*theta(:));
%        tmp = f_r.*f_theta;
%        tmp2(r0<=R) = tmp;
%        tmp3(:, :, lr) = tmp2;
%    end;
%    fn{i+k_max+1} = tmp3;
%end;
