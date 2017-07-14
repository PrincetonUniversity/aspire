% CRYO_PERIODOGRAM Calculate periodogram of images
%
% Usage
%    im_per = cryo_periodogram(im, samples_idx);
%
% Input
%    im: An array of size L1-by-L2-by-n containing images for which we want to
%       calculate the power spectral density.
%    samples_idx: A list of indices in the images to use for the periodogram.
%       If empty, the whole image is used (default empty).
%
% Output
%    im_per: An array of size L1-by-L2-by-n containing the periodograms of the
%       images using the indices in samples_idx. These are in the centered
%       Fourier transform format with the zero frequency at (floor(L1/2)+1,
%       floor(L2/2)+1).

function im_per = cryo_periodogram(im, samples_idx)
    if nargin < 2
        samples_idx = [];
    end

    L = [size(im, 1) size(im, 2)];
    n = size(im, 3);

    if isempty(samples_idx)
        samples_idx = 1:prod(L);
    elseif samples_idx < 1 || samples_idx > prod(L)
        error(['samples_idx must be in the range 1 through ' ...
            'size(im, 1)*size(im 2)']);
    end

    n_samples = numel(samples_idx);

    non_samples_idx = find(~ismember(1:prod(L), samples_idx));

    im = reshape(im, [prod(L) n]);

    im(non_samples_idx,:) = 0;

    im = reshape(im, [L n]);

    im_f = fft2(im);

    im_per = 1/n_samples*abs(im_f).^2;

    im_per = fftshift(fftshift(im_per, 1), 2);
end
