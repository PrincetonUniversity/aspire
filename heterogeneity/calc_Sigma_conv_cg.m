% CALC_SIGMA_CONV_CG Estimate Sigma using conjugate gradient on convolution
%
% Usage
%    [Sigma_est, Sigma_info] = calc_Sigma_conv_cg(kernel_Sigma_f, ...
%        b_Sigma, Sigma_est_opt);
%
% Input
%    kernel_Sigma_f: The Fourier transform of the cubic matrix convolution
%       kernel. This must be a non-centered Fourier transform.
%    b_Sigma: The right-hand side matrix in the normal equations.
%    Sigma_est_opt: An options structure for the estimation of Sigma. Contains
%       the fields:
%          - lambda: A Tikhonov regularization parameter for the least-squares
%             problem. This can be a scalar variable, in which the regular-
%             ization is the identity times lambda. If lambda is a function
%             handle, this is used as a radial filter (default 0).
%          - basis: An N^3-by-M matrix containing as columns vectors the
%             basis in which b_Sigma is described. The resulting output
%             Sigma_est will also be in this basis. Note that the kernel is
%             still in the Euclidean basis.
%          - reweight_matrix: An M-by-M matrix S which reweights the Ln
%             operator so that Ln(Sigma) is replaced by S'*Ln(S*Sigma*S)*S'.
%             Note that this is applied before the basis expansion specified
%             in Sigma_est_opt.basis.
%       This struct is also passed to the conjgrad function for calculating
%       the conjugate gradient solution. See the conjgrad documentation for
%       more details.
%
% Output
%    Sigma_est: The estimated covariance matrix obtained by solving the normal
%       equations with the convolution specified by kernel_Sigma_f and the
%       right-hand side b_Sigma. This is an array of size N-by-N-by-N-by-N-
%       by-N-by-N.
%    Sigma_info: The info structure returned by the conjgrad function. See the
%       conjgrad documentation for more information.

function [Sigma_est, Sigma_info] = calc_Sigma_conv_cg(kernel_Sigma_f, ...
    b_Sigma, Sigma_est_opt)

    if nargin < 3 || isempty(Sigma_est_opt)
        Sigma_est_opt = struct();
    end

    kernel_Sigma_f = volmat_to_vecmat(kernel_Sigma_f);
    b_Sigma = volmat_to_vecmat(b_Sigma);

    Ker_fun = @(Sigma)(apply_conv_Ln(Sigma, kernel_Sigma_f, Sigma_est_opt));

    Ker_fun = @(Sigma)(mat_to_vec(Ker_fun(vec_to_mat(Sigma))));

    [Sigma_est, ~, Sigma_info] = conjgrad(mat_to_vec(b_Sigma), Ker_fun, Sigma_est_opt, []);

    Sigma_est = vec_to_mat(Sigma_est);

    Sigma_est = vecmat_to_volmat(Sigma_est);
end
