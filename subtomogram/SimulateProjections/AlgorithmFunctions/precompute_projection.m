function [Wp, Wm, Crm, nn, fn, r0, InnerP] = precompute_projection(ns, n_r, R, c, L,ss, basis, sample_points)

% Precomputation for simulating projection images for a three-dimensional
% object from coefficients of spherical expansion
%
% INPUT: 
%   ns: number of tilt series, could be 1 for single projection
%   n_r: number of sample points along the radius
%   R: radius of the original volume, (L = 2*R+1)
%   c: compact support radius, usually set to 1/2
%   L: size of the original volume
%   ss: factor for precision requirement, 3 gives machine precision for 
%       Guassian volumes.
%   basis: spherical Bessel basis, output of spherical_expansion.m
%   sample_points: Gaussian quadrature points and weights, from
%       spherical_expansion.m
% OUTPUT:
%   Wp: wigner D matrix for Euler ZYZ angle [-pi/2, pi/2, pi/2], 
%   Wm: wigner D matrix for Euler ZYZ angle [-pi/2, -pi/2, pi/2],
%   Crm: legendre polynomials Ylm(pi/2,0) or with tilt series adjusted to
%       pi/2
%   nn: total number of radius basis used across all frequencies
    %   basis_fb: Fourier Bessel basis with L_max cells, each of size N(l)* 4cR,
    %       as they are sampled on Gaussian quadrature points
%   fn: precomputed inverse function for Fourier Bessel in the real space 
%       and on uniform Cartesian grid
%   r0: radius of the Cartesian grid
%   InnerP: inner products between spherical Bessel basis in 3D and Fourier
%       Bessel basis in 2D

Ctemp = cell(1,L);
Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
Wp = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);
Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
Wm = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);

% precompute legendre polynomials Ylm(pi/2,0), same for all subtomograms
if ns == 1
    theta = pi/2;
else
    tilt = 2/3*pi/(ns-1);
    theta = (-tilt*(ns-1)/2:tilt:tilt*(ns-1)/2)+pi/2;
end
Crm = zeros(L,L,ns); % first index m, second index l
for ll = 0:L-1        
    %%% spherical harmonics for order ll
    m = 0:ll;
    [X, ~, ~]= xlm( ll, m, theta);    
    
    Crm(ll+1,1:ll+1,:) = X;
end
Crm = num2cell(Crm,3);

% precompute fb basis
[ basis_fb, ~] = precomp_fb_norm( n_r, R, c, L,ss ); % checked to third digit

%Computes eigen images, need output from IFT_FB.m.
x1 = -R:1:R;
[ fn, r0 ] = IFT_FBm(R, c, L, x1, ss);

% precompute inner product between 3d and 2d basis
InnerP = FBinnerSB(basis_fb, basis, sample_points); % n of sb by q of fb
nn = length(basis_fb.ang_freqs);

end

