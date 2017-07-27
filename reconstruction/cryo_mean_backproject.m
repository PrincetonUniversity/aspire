% CRYO_MEAN_BACKPROJECT Backproject images for mean volume estimation
%
% Usage
%    im_bp = cryo_mean_backproject(im, params, mean_est_opt);
%
% Input
%    im: An array of size L-by-L-by-n containing images to be backprojected.
%    params: An imaging parameters structure containing the fields:
%          - rot_matrices: A 3-by-3-by-n array containing the rotation
%             matrices of the various projections.
%          - ctf: An L-by-L-by-K array of CTF images (centered Fourier
%             transforms of the point spread functions of the images)
%             containing the CTFs used to generate the images.
%          - ctf_idx: A vector of length n containing the indices in the
%             `ctf` array corresponding to each projection image.
%          - ampl: A vector of length n specifying the amplitude multiplier
%             of each image.
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double'
%             (default) or 'single'.
%
% Output
%    im_bp: The backprojected images, averaged over the whole dataset.

function im_bp = cryo_mean_backproject(im, params, mean_est_opt)
    if nargin < 3 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    L = size(im, 1);
    n = size(im, 3);

    check_imaging_params(params, L, n);

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'double');

    pts_rot = rotated_grids(L, params.rot_matrices);

    im = im_filter(im, params.ctf(:,:,params.ctf_idx));

    im_f = 1/L^2*cfft2(im);

    im_f = bsxfun(@times, im_f, reshape(params.ampl, [1 1 n]));

    if mod(L, 2) == 0
        pts_rot = pts_rot(:,2:end,2:end,:);

        im_f = im_f(2:end,2:end,:);
    end

    pts_rot = reshape(pts_rot, 3, []);
    im_f = reshape(im_f, [], 1);

    im_bp = 1/n*(2/L)*anufft3(im_f, pts_rot, L*ones(1, 3));

    im_bp = real(im_bp);
end
