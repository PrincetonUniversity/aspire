% FB_ROT_COVAR_CTF Rotationally invariant covariance from filtered images
%
% Usage
%    covar_coeff = fb_rot_covar_ctf(coeff, filter_fb, filter_idx, mean_coeff, ...
%       noise_var, basis, covar_est_opt);
%
% Input
%    coeff: The coefficients of images in a 2D Fourier-Bessel basis arranged
%       as a basis.count-by-n array, where n is the number of images.
%    filter_fb: The block diagonal matrices representing filters arranged in a
%       cell array of length K. Each element is a block diagonal matrix
%       represented as another cell array that is manipulated by `blk_diag_*`
%       functions. Each block diagonal matrix corresponds to a basis.count-by-
%       basis-count size matrix.
%    filter_idx: The filter indices of the images in coeff as vector of length
%       n.
%    mean_coeff: The coefficients of the mean image, as obtained from the
%       function fb_rot_mean_ctf.
%    noise_var: The variance of the noise in the images.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%    covar_est_opt: An options structure. Currently, it is sent as options to
%       the conj_grad function. See its documentation for details.
%
% Output
%    covar_coeff: The Fourier-Bessel coefficients of the covariance matrix in
%       the form of cell array representing a block diagonal matrix. These
%       block diagonal matrices may be manipulated using the `blk_diag_*`
%       functions. The covariance is calculated from the images represented
%       by the coeff array, along with all possible rotations and reflections.
%       As a result, the computed covariance matrix is invariant to both
%       reflection and rotation. The effect of the filters in filter_fb are
%       accounted for and inverted to yield a covariance estimate of the
%       unfiltered images.
%
% See also
%    fb_rot_mean_ctf, fb_rot_covar, fb_rot_mean

function covar_coeff = fb_rot_covar_ctf(coeff, ctf_fb, ctf_idx, ...
    mean_coeff, noise_var, basis, covar_est_opt)

    if ~ismember(basis.type, [fb_basis_type() ffb_basis_type()]) || ...
        numel(basis.sz) ~= 2

        error('Basis must be 2D Fourier-Bessel basis.');
    end

    if nargin < 7 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    covar_est_opt = fill_struct(covar_est_opt, ...
        'verbose', 0, ...
        'max_iter', 250, ...
        'rel_tolerance', 1e-12);

    block_partition = blk_diag_partition(ctf_fb{1});

    b = blk_diag_zeros(block_partition, class(coeff));
    A = cell(numel(ctf_fb), 1);
    M = blk_diag_zeros(block_partition, class(ctf_fb{1}{1}));

    for k = unique(ctf_idx(:))'
        coeff_k = coeff(:,ctf_idx == k);

        weight = size(coeff_k, 2)/size(coeff, 2);

        ctf_fb_k = ctf_fb{k};
        ctf_fb_k_t = blk_diag_transpose(ctf_fb_k);

        mean_coeff_k = blk_diag_apply(ctf_fb_k, mean_coeff);

        covar_coeff_k = fb_rot_covar(coeff_k, mean_coeff_k, basis);
        covar_coeff_k = blk_diag_add(covar_coeff_k, ...
            blk_diag_mult(-noise_var, ...
            blk_diag_eye(block_partition, class(covar_coeff_k{1}))));

        b = blk_diag_add(b, ...
            blk_diag_mult(ctf_fb_k_t, ...
            blk_diag_mult(covar_coeff_k, ...
            blk_diag_mult(ctf_fb_k, weight))));

        A{k} = blk_diag_mult(ctf_fb_k_t, ...
            blk_diag_mult(ctf_fb_k, sqrt(weight)));

        M = blk_diag_add(M, A{k});
    end

    cg_opt = covar_est_opt;

    covar_coeff = blk_diag_zeros(block_partition, class(coeff));

    for k = 1:numel(b)
        A_k = cellfun(@(blk)(blk{k}), A, 'uniformoutput', false);
        b_k = b{k};

        S = inv(M{k});

        cg_opt.preconditioner = @(x)(precond(S, x));

        covar_coeff{k} = conj_grad(@(x)(apply(A_k, x)), b_k(:), cg_opt);

        covar_coeff{k} = reshape(covar_coeff{k}, size(A_k{1}, 1)*ones(1, 2));
    end
end

function y = apply(A, x)
    p = size(A{1}, 1);

    x = reshape(x, p*ones(1, 2));

    y = zeros(size(x));

    for k = 1:numel(A)
        y = y + A{k}*x*A{k}';
    end

    y = y(:);
end

function y = precond(S, x)
    p = size(S, 1);

    x = reshape(x, p*ones(1, 2));

    y = S*x*S;

    y = y(:);
end
