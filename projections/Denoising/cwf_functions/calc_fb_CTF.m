function [fb_CTF] = calc_fb_CTF(ctf_rad, Phins_k, sample_points)

% Numerical evaluation of CTF in truncated Fourier Bessel basis
% Phins_k: Radial bessel function, ang freq k
% sample_points: Quadrature points for evaluation, w are the associated
% weights.
% Tejal Bhamre: updates Oct 2015

r=sample_points.r;
w=sample_points.w;
ctf_rad=ctf_rad.*r.*w;

fb_CTF= 2*pi*Phins_k'*(diag(ctf_rad))*Phins_k;
