% CALC_FB_CTF Evaluate radial CTF application in Fourier-Bessel basis
%
% Usage
%    fb_CTF = calc_fb_CTF(ctf_rad, Phins_k, sample_points);
%
% Input
%    ctf_rad: The CTF evaluated at the radial points in sample_points.
%    Phins_k: The radial Fourier-Bessel functions for a fixed angular
%       frequency k. Obtained from precomp_fb.
%    sample_points: The quadrature scheme used. Obtained from precomp_fb.
%
% Output
%    fb_CTF: The matrix corresponding to the application of the CTF in the
%       Fourier-Bessel for the specified radial basis functions.

% Written by Tejal - Oct 2015
% Reformatted and documented by Joakim Andén - 2018-Apr-19

function fb_CTF = calc_fb_CTF(ctf_rad, Phins_k, sample_points)
    r=sample_points.r;
    w=sample_points.w;
    ctf_rad=ctf_rad.*r.*w;

    fb_CTF= 2*pi*Phins_k'*(diag(ctf_rad))*Phins_k;
end
