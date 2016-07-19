function [ basis, sample_points ] = precomp_fb( n_r, R, c )
%%% Description
%The function computes computes the bessel functions
%   Input: 
%          n_r: number of samples on the bessel function at interval [0, c]    
%          c: band limit 
%          R: Compact support radius in real image domain
%   Output:
%          basis.Phi_ns: bessel radial function J(R_{kq}\xi/c) for \xi in [0, c]
%          basis.ang_freqs: the associated angular frequencies (k in the paper)
%          basis.n_theta: the number of samples on concentric rings for polar Fourier transformation.
%          sample_points.r: positions in interval [0, c]
%          sample_points.w: weights
%   Zhizhen Zhao 11/21/2014

%Choose points in [0, c] using Gauss-legendre quadrature rule     
[r, w] = lgwt(n_r, 0, c);
sample_points.r = r;
sample_points.w = w;

%Computes the bessel radial functions
[ basis ] = Bessel_ns_radial(c, R, r);
end

