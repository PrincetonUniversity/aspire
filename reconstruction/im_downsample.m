% IM_DOWNSAMPLE Blur and downsample an image
%
% Usage
%    im_ds = im_downsample(im, N_ds, filter);
%
% Input
%    im: Set of images to be downsampled in the form of an array N-by-N-by-M,
%       where M is the number of images.
%    N_ds: The desired resolution of the downsampled images.
%    filter: The type of filter to use: 'gaussian', or 'sinc' (default
%       'gaussian').
%
% Output
%    im_ds: An array of the form N_ds-by-N_ds-by-M consisting of the blurred
%       and downsampled images.

function im_ds = im_downsample(im, N_ds, filter)
    if nargin < 3 || isempty(filter)
        filter = 'gaussian';
    end

    N = size(im, 1);

    mesh = mesh_2d(N);
    mesh_ds = mesh_2d(N_ds);

    im_ds = zeros([N_ds*ones(1, 2) size(im, 3)]);

    % NOTE: Cannot use interp2 because it's not compatible with ndgrid.
    ds_fun = @(im)(interpn(mesh.x, mesh.y, im, mesh_ds.x, mesh_ds.y, ...
        'linear', 0));

    if strcmp(filter, 'gaussian')
        c = N/N_ds/(2*sqrt(2*log(2)));
        blur = exp(-[-ceil(3*c):ceil(3*c)].^2/(2*c^2));

        for s = 1:size(im_ds, 3)
            im_s = conv2(blur, blur, im(:,:,s), 'same');
            im_ds(:,:,s) = ds_fun(im_s);
        end
    else
        mask = (abs(mesh.x) < N_ds/N) & (abs(mesh.y) < N_ds/N);

        for s = 1:size(im_ds, 3)
            im_ds(:,:,s) = ds_fun(im(:,:,s));
        end
    end
end
