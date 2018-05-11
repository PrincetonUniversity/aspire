% CRYO_CONJ_GRAD_MEAN Solve for mean volume using conjugate gradient
%
% Usage
%    mean_est = cryo_conj_grad_mean(kernel_f, im_bp, basis, ...
%       precond_kernel_f, mean_est_opt);
%
% Input
%    kernel_f: The centered Fourier transform of the projection-backprojection
%       operator obtained from `cryo_mean_kernel_f`.
%    im_bp: An array of size L-by-L-by-L containing the backprojected images
%       obtained from `cryo_mean_backproject`.
%    basis: A basis object used for representing the volumes, such as
%       dirac_basis(L*ones(1, 3)).
%    precond_kernel_f: If present, the Fourier transform of a kernel that is
%       used to precondition the projection-backprojection operator (default
%       empty).
%    mean_est_opt: An options structure. No options are used in the function
%       itself, but it is passed on to the `conj_grad` function, so any options
%       for that function should be specified here.
%
% Output
%    mean_est: An array of size L-by-L-by-L containing the least-squares
%       estimate obtained by solving the equation A*x = b, where A is the linear
%       mapping represented by convolving with `kernel_f` and b are the
%       backprojected images `im_bp`. The equation is solved using the
%       conjugate gradient method.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

function [mean_est, cg_info] = cryo_conj_grad_mean(kernel_f, im_bp, basis, ...
    precond_kernel_f, mean_est_opt)

    if nargin < 4
        precond_kernel_f = [];
    end

    if nargin < 5 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    L = size(im_bp, 1);

    if ndims(im_bp) ~= 3 || size(im_bp, 2) ~= L || size(im_bp, 3) ~= L
        error('Input `im_bp` must be an array of size L-by-L-by-L.');
    end

    if ~is_basis(basis) || any(basis.sz ~= L*ones(1, 3))
        error(['Input `basis` must be a basis object representing ' ...
            'volumes of size L-by-L-by-L.']);
    end

    fun = @(vol_basis)( ...
        apply_mean_kernel(vol_basis, kernel_f, basis, mean_est_opt));

    if ~isempty(precond_kernel_f)
        precond_fun = @(vol_basis)( ...
            apply_mean_kernel(vol_basis, precond_kernel_f, basis, ...
            mean_est_opt));

        mean_est_opt.preconditioner = precond_fun;
    end

    im_bp_basis = basis.evaluate_t(im_bp);

    [mean_est_basis, ~, cg_info] = conj_grad(fun, im_bp_basis, mean_est_opt);

    mean_est = basis.evaluate(mean_est_basis);
end
