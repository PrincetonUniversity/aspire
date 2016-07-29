% VOL_PROJECT Project volume along rotation
%
% Usage
%    im = vol_project(vol, rot_matrices);
%
% Input
%    vol: An N-by-N-by-N-by-K array containing the voxel structures of K
%       volumes. If K = 1, the volumes are broadcast over the set of projection
%       directions.
%    rot_matrices: A set of rotation matrices of the form 3-by-3-by-K, corre-
%       sponding to K different projection directions. If K = 1, the projec-
%       tion directions are broadcast over the set of volumes.
%    full_proj: Specifies whether to do a full projection or restrict the to
%       the unit ball in the Fourier domain (default false).
%
% Output
%    im: An N-by-N-by-K array containing the projections of the volumes in the
%       specified directions.

function im = vol_project(vol, rot_matrices, full_proj)
    global g_parallel;

    if nargin < 3 || isempty(full_proj)
        full_proj = false;
    end

    N = size(vol, 1);

    if size(rot_matrices,2) == 1
        error('second argument must be a set of rotation matrices');
    end

    if ndims(vol) == 2
        warning('calling with vol vectors is deprecated');
        im = vol_project(vec_to_vol(vol), rot_matrices);
        return;
    end

    n_vols = size(vol, 4);
    n_rots = size(rot_matrices, 3);

    n = max(n_vols, n_rots);

    if (n_vols ~= 1 && n_vols ~= n) || (n_rots ~= 1 && n_rots ~= n)
        error(['number of volumes and rotations need to be either the ' ...
            'same or one']);
    end

    vol_ind = broadcast_idx(n_vols, n);
    rot_ind = broadcast_idx(n_rots, n);

    if ~full_proj
        [pts_rot, mask] = rotated_grids(N, rot_matrices);
    else
        [pts_rot, mask] = rotated_grids(N+1, rot_matrices, false, true);
    end
    n_pts = sum(mask(:));

    % NOTE: Cannot define it as N^2 by n because MATLAB complains in parfor 
    % otherwise (why?)
    im_f_masked = zeros([n_pts n], class(vol));

    if ~isempty(g_parallel) && g_parallel
        parfor s = 1:n
            im_f_masked(:,s) = 2/N*nufft3(vol(:,:,:,vol_ind(s)), ...
                pi*pts_rot(:,:,rot_ind(s))');
        end
    else
        for s = 1:n
            im_f_masked(:,s) = 2/N*nufft3(vol(:,:,:,vol_ind(s)), ...
                pi*pts_rot(:,:,rot_ind(s))');
        end
    end

    if ~full_proj
        im_f = zeros([N^2 n], class(im_f_masked));
        im_f(mask,:) = im_f_masked;
        im_f = vec_to_im(im_f);
    else
        im_f = reshape(im_f_masked, [(N+1)*ones(1, 2) n]);
        im_f(1,:,:) = (im_f(1,:,:)+im_f(N+1,:,:))/2;
        im_f(:,1,:) = (im_f(:,1,:)+im_f(:,N+1,:))/2;
        im_f = im_f(1:N,1:N,:);
    end

    im = real(centered_ifft2(im_f));
end
