% FBCOEFF_NFFT Compute Fourier-Bessel expansion of images
%
% Usage
%    coeff_pos_k = FBcoeff_nfft(data, R, basis, sample_points)
%
% Input
%    data: A cell array containing the images. Each element is of size
%       L-by-L-by-n, where each image is L-by-L and there are n of them.
%    R: The radius of the disk on which the basis is supported.
%    basis: The basis struct, obtained from precomp_fb.
%    sample_points: The struct describing the quadrature, obtained from
%        precomp_fb.
%
% Output
%    coeff_pos_k: A cell array of the Fourier-Bessel coefficients. Each cell
%       contains the coefficients for a fixed angular frequency, ranging from
%       0 to max(basis.ang_freqs)). These coefficients are two-dimensional
%       arrays, with the first dimension corresponding to radial frequency and
%       the second corresponding to image index. The images are numbered
%       as though data were concatenated along the third dimension.
%
% Description
%    The coefficients are obtained using a polar Fourier transform on the
%    image, followed by a one-dimensional FFT along the angular direction
%    and an inner product with the radial Fourier-Bessel basis functions.
%    In the notation of Bessel_ns_radial, the coefficient corresponding
%    to the basis vector \phi_{k, q}(r) e^{i k \theta) is found in the
%    array coeff_pos_k{k}(q,:).

% Written by Zhizhen Zhao - 09/2015
% Reformatted and documented by Joakim Anden - 2018-Apr-13

function coeff_pos_k = FBcoeff_nfft(data, R, basis, sample_points)
    % Initial image size
    L = size(data{1}, 1);
    orig = floor(L/2)+1;
    num_pool = numel(data);

    % Read input data
    Phi_ns = basis.Phi_ns;
    ang_freqs = basis.ang_freqs;
    n_theta = basis.n_theta;
    clear basis;

    max_ang_freqs = max(ang_freqs);

    w = sample_points.w;
    r = sample_points.r;
    w = r.*w;

    freqs = pft_freqs(r, n_theta);
    Precomp.n_theta = n_theta;
    Precomp.n_r = length(r);
    Precomp.freqs = freqs;

    scale = 2*pi/(n_theta);

    % Evaluate expansion coefficients
    coeff_pos_k = cell(max_ang_freqs+1, 1);
    pos_k = cell(num_pool, 1);

    parfor i = 1:num_pool
        tmp = data{i};
        tmp = tmp(orig-R:orig+R-1,orig-R:orig+R-1,:);
        tmp2 = cryo_pft_nfft(tmp, Precomp);

        % 1D FFT on concentric rings
        pf_f = scale*fft(tmp2, [], 2);
        pos_k{i} = pf_f(:,1:max(ang_freqs)+1,:);
    end

    pos_k = cat(3, pos_k{:});

    clear pf_f tmp tmp2;

    for i = 1:max_ang_freqs+1
        tmp = squeeze(pos_k(:,i,:));
        coeff_pos_k{i} = Phi_ns{i}.'*diag(w)*tmp;
    end
end
