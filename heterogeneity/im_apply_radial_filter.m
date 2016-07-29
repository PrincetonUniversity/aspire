% IM_APPLY_RADIAL_FILTER Mask image in Fourier and apply radial filter
%
% Usage
%    im = im_apply_radial_filter(im, filter_fun);
%
% Input
%    im: A collection of images in an array of size N-by-N-by-n, where N is
%       the resolution of the images and n is the number of images.
%    filter_fun: A single function handle whose input is a set of radial
%       frequencies and whose output is the value of the filter transfer
%       function at those radii. Alternatively, filter_fun can be a cell array
%       of function handles. In this case, the array must be of length n and
%       each image is then filtered by the corresponding filter function.
%
% Output
%    im: The images filtered by the single filter given by filter_fun or by
%       the corresponding filter in the case of multiple function handles.

function im = im_apply_radial_filter(im, filter_fun)
    N = size(im,1);

    [im, sz_roll] = unroll_dim(im, 3);

    mesh2d = mesh_2d(N);
    im_f = centered_fft2(im);
    im_f = im_to_vec(im_f);
    im_f = bsxfun(@times, im_f, mesh2d.mask(:));

    if iscell(filter_fun) && numel(filter_fun) == 1
        filter_fun = filter_fun{1};
    end

    [r_vals,~,ind] = unique(mesh2d.r(:));
    if ~iscell(filter_fun)
        f_vals = filter_fun(r_vals);
        im_f(mesh2d.mask,:) = bsxfun(@times, im_f(mesh2d.mask,:), f_vals(ind));
    else
        if length(filter_fun) ~= size(im, 3) && length(filter_fun) ~= 1
            error('number of filters must match number of images or equal 1');
        end

        for s = 1:size(im, 3)
            f_vals = filter_fun{s}(r_vals);
            im_f(mesh2d.mask,s) = im_f(mesh2d.mask,s).*f_vals(ind);
        end
    end

    im_f = vec_to_im(im_f);
    im = centered_ifft2(im_f);

    im = roll_dim(im, sz_roll);
end
