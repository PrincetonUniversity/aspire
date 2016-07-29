% CALC_MU_CONV_CG Estimate mu using conjugate gradient on convolution
%
% Usage
%    [mu_est, mu_info] = calc_mu_conv_cg(kernel_mu_f, b_mu, mu_est_opt);
%
% Input
%    kernel_mu_f: The Fourier transform of the cubic convolution kernel. This
%       must be a non-centered Fourier transform.
%    b_mu: The right-hand side volume in the normal equations.
%    mu_est_opt: An options structure for the estimation of mu. Contains
%       the fields:
%          - lambda: A Tikhonov regularization parameter for the least-squares
%             problem. This can be a scalar variable, in which the regular-
%             ization is the identity times lambda. If lambda is a function
%             handle, this is used as a radial filter (default 0).
%          - basis: An N^3-by-M matrix containing as columns vectors the
%             basis in which b_mu is described. The resulting output mu_est
%             will also be in this basis. Note that the kernel is still in
%             the Euclidean basis.
%       This struct is also passed to the conjgrad function for calculating
%       the conjugate gradient solution. See the conjgrad documentation for
%       more details.
%
% Output
%    mu_est: The estimated volume obtained by solving the normal equations with
%       the convolution specified by kernel_Sigma_f and the right-hand side
%       b_Sigma.
%    mu_info: The info structure returned by the conjgrad function. See the
%       conjgrad documentation for more information.

function [mu_est, mu_info] = calc_mu_conv_cg(kernel_mu_f, b_mu, mu_est_opt)
    if nargin < 3 || isempty(mu_est_opt)
        mu_est_opt = struct();
    end

    mu_est_opt = fill_struct(mu_est_opt, ...
        'lambda', 0, ...
        'basis', []);

    if ~isempty(mu_est_opt.basis)
        N = round(size(mu_est_opt.basis, 1)^(1/3));
    else
        N = round(size(b_mu, 1)^(1/3));
    end

    vol_fun = @(x)(x);
    if ~isempty(mu_est_opt.basis)
        vol_fun = @(x)(mu_est_opt.basis*x);
    end

    Ker_fun = vol_fun;

    Ker_fun = @(x)(vol_to_vec(conv_vol(vec_to_vol(Ker_fun(x)), kernel_mu_f)));

    if isnumeric(mu_est_opt.lambda)
        Ker_fun = @(x)(Ker_fun(x) + mu_est_opt.lambda*vol_fun(x));
    elseif isa(mu_est_opt.lambda, 'function_handle');
        reg_fun = @(x)(vol_to_vec(vol_apply_radial_filter( ...
            vec_to_vol(vol_fun(x)), mu_est_opt.lambda)));

        Ker_fun = @(x)(Ker_fun(x) + reg_fun(x));
    end

    if ~isempty(mu_est_opt.basis)
        Ker_fun = @(x)(mu_est_opt.basis'*Ker_fun(x));
    end

    b_mu = vol_to_vec(b_mu);

    [mu_est, ~, mu_info] = conjgrad(b_mu, Ker_fun, mu_est_opt);

    mu_est = real(vec_to_vol(mu_est));
end

