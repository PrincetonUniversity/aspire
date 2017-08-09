% IM_DOWNSAMPLE Blur and downsample an image
%
% Usage
%    im_ds = im_downsample(im, N_ds);
%
% Input
%    im: Set of images to be downsampled in the form of an array N-by-N-by-M,
%       where M is the number of images.
%    N_ds: The desired resolution of the downsampled images.
%
% Output
%    im_ds: An array of the form N_ds-by-N_ds-by-M consisting of the blurred
%       and downsampled images.

function im_ds = im_downsample(im, N_ds)
    N = size(im, 1);

    c = N/N_ds/(2*sqrt(2*log(2)));
    blur = exp(-[-ceil(3*c):ceil(3*c)].^2/(2*c^2));

    mesh = mesh_2d(N);
    mesh_ds = mesh_2d(N_ds);

    im_ds = zeros([N_ds*ones(1, 2) size(im, 3)]);

    for s = 1:size(im_ds, 3)
        im_s = conv2(blur, blur, im(:,:,s), 'same');
        % NOTE: Cannot use interp2 because it's not compatible with ndgrid.
        im_ds(:,:,s) = interpn(mesh.x, mesh.y, im_s, mesh_ds.x, mesh_ds.y, ...
            'linear', 0);
    end
end
