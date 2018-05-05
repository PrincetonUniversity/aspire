% IM_FILTER Filter image by a transfer function
%
% Usage
%    im = im_filter(im, filter_f);
%
% Input
%    im: An array of size L-by-L-by-n containing images to be filtered.
%    filter_f: The centered Fourier transform of a filter (its transfer
%       function) or a set of filters. The first two dimensions of this
%       filter must be equal to the first two dimensions of `im`. If one
%       filter is given, it is applied to all images. If not, the shape of
%       `filter_f` must be compatible with `im` and will be used to filter
%       each image individually.
%
% Output
%    im_filtered: The filtered images.

function im_filtered = im_filter(im, filter_f)
    n_im = size(im, 3);

    n_filter = size(filter_f, 3);

    if n_filter ~= 1 && n_filter ~= n_im
        error(['The number of filters must either be 1 or match the ' ...
            'number of images.']);
    end

    if size(im, 1) ~= size(filter_f, 1) || ...
        size(im, 2) ~= size(filter_f, 2)
        error('The size of the images and filters must match.');
    end

    im_f = cfft2(im);

    im_filtered_f = bsxfun(@times, im_f, filter_f);

    im_filtered = icfft2(im_filtered_f);
end
