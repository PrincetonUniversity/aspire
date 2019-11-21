function [ basis ]=Bessel_ns_radial_norm_par( c, R, r, L, ss )
%%% Description
%The function computes the Fourier-Bessel Basis with positive angular frequencies
%   This evaluate Bessel radial functions for an image of size
%   Table bessel.mat gives the R_ns (3rd column) for n (1st column) and s (2nd column)
% Input:  
%       1. c: band limit
%       2. R: compact support radius in real domain
%	    3. r: poistions on [0, 1]
% Output: 
%       1. basis.Phi_ns: Bessel radial functions in cell structure
%       2. basis.ang_freqs: angular frequencies
%       3. basis.rad_freqs: radial frequencies
%       4. basis.n_theta number of samples on concentric ring          
% Zhizhen Zhao 11/25/2014
% Yuan Liu 8/9/2017 modify the normalization

% cell i in basis.Phi_ns stands for angular frequency i-1, in each matrix
% inside cell each column stands for a radial frequency

load bessel.mat  %table for Bessel zeros that   
B = bessel(bessel(:, 4)<= ss*pi*c*R & bessel(:,1)<=L, :);  %Choose functions that R_{k, q+1} \leq \pi N
clear bessel
ang_freqs = B(:, 1);
max_ang_freqs = max(ang_freqs);
%n_theta = ceil(16*c*R);
n_theta = max_ang_freqs+1;
if mod(n_theta, 2)==0
   n_theta = n_theta +1;
end;
ang_freqs = B(:, 1);
rad_freqs = B(:, 2);
R_ns = B(:, 3);
Phi_ns=zeros(length(r), size(B, 1));
Phi = cell(max_ang_freqs+1, 1);

%For matlab2014 above, use parpool() and delete(gcp)
%parpool();
%matlabpool open;
parfor i=1:size(B, 1)
    r0=r*R_ns(i)/c;
    [ F ]=besselj(ang_freqs(i), r0); %Bessel radial functions
    tmp = besselj(ang_freqs(i)+1, R_ns(i))^2;
    %pi*besselj(ang_freqs(i)+1, R_ns(i))^2; %pi*(pi*besselj(ang_freqs(i)+1.5, R_ns(i))^2)./(2*R_ns(i)); change for 3D
    %normalization that \int_0^1 Phi_ns(r) Phi_ns(r) r dr = pi/2
    Phi_ns(:, i)=1/(c*sqrt(tmp/2))*F; 
end;

%put it in the cell structure
parfor i=1:max_ang_freqs+1
    Phi{i} = Phi_ns(:, ang_freqs == i-1);
end;
%delete(gcp)
%matlabpool close;

basis.Phi_ns = Phi;
basis.ang_freqs = ang_freqs;
basis.rad_freqs = rad_freqs;
basis.n_theta = n_theta;





end
 