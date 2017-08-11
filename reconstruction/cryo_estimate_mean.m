% CRYO_ESTIMATE_MEAN Estimate mean using least squares and conjugate gradients
%
% Usage
%    mean_est = cryo_estimate_mean(im, params, basis, mean_est_opt);
%
% Input
%    im: An array of size L-by-L-by-n containing projection images.
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
%    basis: A basis object used for representing the volumes (default
%       dirac_basis(L*ones(1, 3))).
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double'
%             (default) or 'single'.
%          - 'preconditioner': One of the following values specifying the
%             preconditioner for the conjugate gradient method:
%                - 'none': No preconditioner is used.
%                - 'circulant': Uses `circularize_kernel_f` to obtain a
%                   circulant approximation to the projection-backprojection
%                   kernel whose inverse is then used as a preconditioner
%                   (default).
%                - function handle: In this case, the function handle is passed
%                   on directly to the `conj_grad` function.
%       The struct is also passed on to the `conj_grad` function, so any options
%       to that function should be passed here.
%
% Output
%    mean_est: The estimated mean volume, in the form of an L-by-L-by-L array.
%       It minimizes the objective
%
%          1/n sum_{s=1}^n |P_s x - y_s|^2
%
%       where x is the volume, P_s are the imaging mappings (projection,
%       CTF filtering, multiplication), and y_s are the observed images. This
%       is achieved by forming the normal equations and solving them using the
%       conjugate gradient method.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

function [mean_est, cg_info] = cryo_estimate_mean(im, params, basis, ...
    mean_est_opt)

    if nargin < 3
        basis = [];
    end

    if nargin < 4 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    L = size(im, 1);
    n = size(im, 3);

    check_imaging_params(params, L, n);

    mean_est_opt = fill_struct(mean_est_opt, ...
        'preconditioner', 'circulant', ...
        'precision', 'double');

    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kernel_f = cryo_mean_kernel_f(L, params, mean_est_opt);

    precond_kernel_f = [];

    if ischar(mean_est_opt.preconditioner)
        if strcmp(mean_est_opt.preconditioner, 'none')
            precond_kernel_f = [];
        elseif strcmp(mean_est_opt.preconditioner, 'circulant')
            precond_kernel_f = 1./circularize_kernel_f(kernel_f);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        mean_est_opt.preconditioner = @(x)(x);
    end

    im_bp = cryo_mean_backproject(im, params, mean_est_opt);

    [mean_est, cg_info] = cryo_conj_grad_mean(kernel_f, im_bp, basis, ...
        precond_kernel_f, mean_est_opt);
end
