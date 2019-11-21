function [ basis, maxL ]=Bessel_ns_radial_sphere( c, R, r, ss, maxL )
%%% Description
% The function computes the Fourier-Bessel Basis with positive angular frequencies
%   This evaluate Bessel radial functions for an image of size
%   Table bessel.mat gives the R_ns (3rd column) for n (1st column) and s (2nd column)
% Input:  
%   c: band limit
%   R: compact support radius in real domain
%	r: poistions on [0, 1]
%   L: size of the original volume
%   ss: factor for precision requirement, 3 gives machine precision for 
%       Guassian volumes.
% Output: 
%   basis.Phi_ns: Bessel radial functions in cell structure
%   basis.ang_freqs: angular frequencies
%   basis.rad_freqs: radial frequencies
%   basis.n_theta: number of samples on concentric ring          
% Zhizhen Zhao 11/25/2014
% Yuan Liu 8/9/2017 modify the normalization

% cell i in basis.Phi_ns stands for angular frequency i-1, in each matrix
% inside cell each column stands for a radial frequency

load SphericalBesselL700R300pi.mat  %table for Bessel zeros that  

if ~exist('maxL','var')
    B = bessel(bessel(:, 4)<= ss*pi*c*R, :); % c*R
    maxL = max(B(:,1)) + 1;
else
    B = bessel(bessel(:, 4)<= ss*pi*c*R & bessel(:,1)<maxL, :);  %Choose functions that R_{k, q+1} \leq 2*pi*c*R
end
clear bessel
ang_freqs = B(:, 1);
max_ang_freqs = max(ang_freqs);

ang_freqs = B(:, 1);
rad_freqs = B(:, 2);
R_ns = B(:, 3);
Phi_ns=zeros(length(r), size(B, 1));
Phi = cell(max_ang_freqs+1, 1);

for i=1:size(B, 1) %par
    r0=r*R_ns(i)/(c); %k_nl*r
    [ F ]=besselj(ang_freqs(i)+1/2, r0).*sqrt(pi./(2*r0)); % spherical Bessel radial functions
    tmp = (besselj(ang_freqs(i)+3/2, R_ns(i))*sqrt(pi/(2*R_ns(i))))^2;
    Phi_ns(:, i)=F/sqrt((c)^3/2*tmp); 
end

%put it in the cell structure
for i=1:max_ang_freqs+1 %par
    Phi{i} = Phi_ns(:, ang_freqs == i-1);
end;
%delete(gcp)

basis.Phi_ns = Phi;
basis.ang_freqs = ang_freqs;
basis.rad_freqs = rad_freqs;

end
 
