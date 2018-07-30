% PRECOMP_FB Precompute Fourier-Bessel basis
%
% Usage
%    [basis, sample_points] = precomp_fb(n_r, R, c);
%
% Input
%    n_r: Number of samples on the interval [0, c].
%    R: Radius of the disk on which the basis is supported.
%    c: Bandlimit in frequency.
%
% Output
%    basis: A struct containing the fields
%          - Phi_ns: cell array of Bessel radial functions,
%          - ang_freqs: angular frequencies,
%          - rad_freqs: radial frequencies, and
%          - n_theta: number of samples on concentric rings.
%    sample_points: A struct containing the fields
%          - r: quadrature nodes on [0, c],
%          - w: the corresponding quadrature weights.
%
% Description
%    The function computes the Fourier-Bessel basis for a given number of
%    nodes (n_r), support size (R), and bandlimit (c). The radial basis
%    functions are sampled on the n_r nodes sample_points.r between 0 and c
%    in frequency. For more information on the basis structure, see the
%    Bessel_ns_radial documentation.

% Written by Zhizhen Zhao - 11/21/2014
% Reformatted and documented by Joakim Anden - 2018-Apr-13

function [basis, sample_points] = precomp_fb(n_r, R, c)
    % Choose points in [0, c] using Gauss-legendre quadrature rule     
    [r, w] = lgwt(n_r, 0, c);
    sample_points.r = r;
    sample_points.w = w;

    % Computes the bessel radial functions
    basis = Bessel_ns_radial(c, R, r);
end
