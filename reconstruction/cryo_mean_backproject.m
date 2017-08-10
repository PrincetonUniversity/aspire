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
%          - shifts: An array of size 2-by-n containing the offsets in x and y
%             y (the second and first dimensions of im) which were applied
%             after the projection.
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double'
%             (default) or 'single'.
%          - 'batch_size': The size of the batches in which to compute the
%             backprojection, if set to a non-empty value. If empty, there are
%             no batchces and the entire set of images is used. A small batch
%             size can help with NUFFT libraries (like the Chemnitz NFFT can
%             help with certain package), which cannot handle too many nodes
%             at once (default empty).
%
% Output
%    im_bp: The backprojected images, averaged over the whole dataset.

function im_bp = cryo_mean_backproject(im, params, mean_est_opt)
    if nargin < 3 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    if size(im, 1) ~= size(im, 2) || size(im, 1) == 1 || ndims(im) > 3
        error('Input `im` must be of size L-by-L-by-n for L > 1.');
    end

    L = size(im, 1);
    n = size(im, 3);

    check_imaging_params(params, L, n);

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'double', ...
        'batch_size', []);

    if ~isempty(mean_est_opt.batch_size)
        % To do batch, simply take a subset of the parameters, call the
        % function for those subsets, and sum.
        batch_size = mean_est_opt.batch_size;

        mean_est_opt.batch_size = [];

        batch_ct = ceil(n/batch_size);

        im_bp = zeros(L*ones(1, 3), mean_est_opt.precision);

        for batch = 1:batch_ct
            s1 = (batch-1)*batch_size+1;
            s2 = min(batch*batch_size, n);

            batch_params = subset_params(params, s1:s2);
            batch_im = im(:,:,s1:s2);

            batch_im_bp = cryo_mean_backproject(batch_im, batch_params, ...
                mean_est_opt);

            im_bp = im_bp + (s2-s1+1)/n*batch_im_bp;
        end

        return;
    end

    pts_rot = rotated_grids(L, params.rot_matrices);

    im = im_translate(im, -params.shifts);

    im = im_filter(im, params.ctf(:,:,params.ctf_idx));

    im = permute(im, [2 1 3]);

    im_f = 1/L^2*cfft2(im);

    im_f = bsxfun(@times, im_f, reshape(params.ampl, [1 1 n]));

    if mod(L, 2) == 0
        pts_rot = pts_rot(:,2:end,2:end,:);

        im_f = im_f(2:end,2:end,:);
    end

    pts_rot = reshape(pts_rot, 3, []);
    im_f = reshape(im_f, [], 1);

    im_bp = 1/n*anufft3(im_f, pts_rot, L*ones(1, 3));

    im_bp = real(im_bp);
end
