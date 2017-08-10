% IM_BACKPROJECT Backproject images along rotation
%
% Usage
%    vol = im_backproject(im, rot_matrices)
%
% Input
%    im: An L-by-L-by-n array of images to backproject.
%    rot_matrices: An 3-by-3-by-n array of rotation matrices corresponding to
%       viewing directions.
%
% Output
%    vol: An L-by-L-by-L-by-n array of volumes corresponding to the desired
%       backprojections.

function vol = im_backproject(im, rot_matrices)
    L = size(im, 1);

    if size(im, 2) ~= L
        error('im must be of the form L-by-L-by-K');
    end

    n = size(im, 3);

    if size(rot_matrices, 3) ~= n
        error(['The number of rotation matrices must match the number ' ...
            'of images.']);
    end

    pts_rot = rotated_grids(L, rot_matrices);

    pts_rot = reshape(pts_rot, [3 L^2*n]);

    im = permute(im, [2 1 3]);

    im_f = 1/L^2*cfft2(im);

    if mod(L, 2) == 0
        im_f(1,:,:) = 0;
        im_f(:,1,:) = 0;
    end

    im_f = reshape(im_f, [L^2*n 1]);

    vol = anufft3(im_f, pts_rot, L*ones(1, 3));

    vol = real(vol);
end
