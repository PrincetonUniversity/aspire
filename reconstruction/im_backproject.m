% IM_BACKPROJECT Backproject images along rotation
%
% Usage
%    vol = im_backproject(im, rot_matrices, half_pixel);
%
% Input
%    im: An L-by-L-by-n array of images to backproject.
%    rot_matrices: An 3-by-3-by-n array of rotation matrices corresponding to
%       viewing directions.
%    half_pixel: If true, centers the rotation around a half-pixel (default false).
%
% Output
%    vol: An L-by-L-by-L-by-n array of volumes corresponding to the desired
%       backprojections.

function vol = im_backproject(im, rot_matrices, half_pixel)
    if nargin < 3 || isempty(half_pixel)
        half_pixel = false;
    end

    L = size(im, 1);

    if size(im, 2) ~= L
        error('im must be of the form L-by-L-by-K');
    end

    n = size(im, 3);

    if size(rot_matrices, 3) ~= n
        error(['The number of rotation matrices must match the number ' ...
            'of images.']);
    end

    if mod(L, 2) == 0 && half_pixel
        [X, Y] = ndgrid(-L/2:L/2-1, -L/2:L/2-1);
    end

    pts_rot = rotated_grids(L, rot_matrices, half_pixel);

    pts_rot = reshape(pts_rot, [3 L^2*n]);

    im = permute(im, [2 1 3]);

    if mod(L, 2) == 0 && half_pixel
        phase_shift = 2*pi*(X+Y)/(2*L);
        im = bsxfun(@times, im, exp(-1i*phase_shift));
    end

    im_f = 1/L^2*cfft2(im);

    if mod(L, 2) == 0
        if ~half_pixel
            im_f(1,:,:) = 0;
            im_f(:,1,:) = 0;
        else
            phase_shift = -reshape(sum(pts_rot, 1), [L*ones(1, 2) n])/2;
            phase_shift = phase_shift + 2*pi*(X+Y+1)/(2*L);
            im_f = bsxfun(@times, im_f, exp(-1i*phase_shift));
        end
    end

    im_f = reshape(im_f, [L^2*n 1]);

    vol = anufft3(im_f, pts_rot, L*ones(1, 3));

    vol = real(vol);
end
