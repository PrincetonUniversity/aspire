% BESSEL_NS_RADIAL Compute Fourier-Bessel basis, positive angular frequencies
%
% Usage
%    basis = Bessel_ns_radial(c, R, r);
%
% Input
%    c: The band limit in frequency.
%    R: The radius of the disk on which the basis is supported.
%    r: Radial sample points in frequency between 0 and c.
%
% Output
%    basis: A struct containing the fields:
%          - Phi_ns: cell array of Bessel radial functions,
%          - ang_freqs: angular frequencies,
%          - rad_freqs: radial frequencies, and
%          - n_theta: number of samples on concentric rings.
%
% Description
%    The radial basis functions computed are of the form
%
%       \phi_{k, q}(u) = c_{k, q} J_k ( R_{k, q} u / c ),
%
%    where u is the radial frequency in [0, c], J_k is the kth Bessel
%    function, R_{k, q} is its qth zero, and c_{k, q} is a normalization
%    constant so that
%
%       \int_0^c | \phi_{k, q} (u) |^2 u du = 1/(2 pi).
%
%    Its values are sampled at the points given by r.
%
%    The basis.Phi_ns cell array contains these sampled radial basis functions
%    arranged according to the angular frequency k. Each element is a matrix
%    whose columns correspond to the radial frequency q.

% Written by Zhizhen Zhao - 11/25/2014
% Reformatted and documented by Joakim Anden - 2018-Apr-13

function basis = Bessel_ns_radial(c, R, r)

    % Choose functions that R_{k, q+1} \leq \pi N
    f = load(fullfile(aspire_root(), 'projections', 'Denoising', 'ffb', ...
        'bessel.mat'));
    B = f.bessel(f.bessel(:, 4)<= 2*pi*c*R, :);
    clear f;

    % First column of B contains Bessel order k
    % Second column contains zero index q
    % Third column contains location of qth zero R_{k, q} of J_k(u)
    % Fourth column contains R_{k, q+1}

    ang_freqs = B(:, 1);
    max_ang_freqs = max(ang_freqs);
    n_theta = ceil(16*c*R);
    if mod(n_theta, 2)==1
       n_theta = n_theta +1;
    end
    ang_freqs = B(:, 1);
    rad_freqs = B(:, 2);
    R_ns = B(:, 3);
    Phi_ns=zeros(length(r), size(B, 1));
    Phi = cell(max_ang_freqs+1, 1);

    parfor i=1:size(B, 1)
        r0=r*R_ns(i)/c;
        F=besselj(ang_freqs(i), r0);
        tmp = pi*besselj(ang_freqs(i)+1, R_ns(i))^2;
        % Normalization \int_0^1 Phi_ns(r) Phi_ns(r) r dr = pi/2
        Phi_ns(:, i)=1/(c*sqrt(tmp))*F;
    end

    parfor i=1:max_ang_freqs+1
        Phi{i} = Phi_ns(:, ang_freqs == i-1);
    end

    basis.Phi_ns = Phi;
    basis.ang_freqs = ang_freqs;
    basis.rad_freqs = rad_freqs;
    basis.n_theta = n_theta;
end
