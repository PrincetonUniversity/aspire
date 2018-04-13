% BESSEL_NS_RADIAL Compute Fourier-Bessel basis, positive angular frequencies
%
% Usage
%    basis = Bessel_nr_radial(c, R, r);
%
% Input
%    c: The band limit in frequency.
%    R: The radius of the disk on which the basis is supported.
%    r: Radial sample points in frequency between 0 and c.
%
% Output
%    basis: A struct containing the fields:
%           - Phi_ns: cell array of Bessel radial functions,
%           - ang_freqs: angular frequencies,
%           - rad_freqs: radial frequencies, and
%           - n_theta: number of samples on concentric rings.

% Written by Zhizhen Zhao - 11/25/2014
% Reformated and documented by Joakim Anden - 2018-Apr-13

function basis = Bessel_ns_radial(c, R, r)

    % Choose functions that R_{k, q+1} \leq \pi N
    f = load(fullfile(aspire_root(), 'projections', 'Denoising', 'ffb', ...
        'bessel.mat'));
    B = f.bessel(f.bessel(:, 4)<= 2*pi*c*R, :);
    clear f;

    ang_freqs = B(:, 1);
    max_ang_freqs = max(ang_freqs);
    n_theta = ceil(16*c*R);
    if mod(n_theta, 2)==0
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
