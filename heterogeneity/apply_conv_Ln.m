% APPLY_CONV_LN Apply Ln using convolution
%
% Usage
%    volmat = apply_conv_Ln(volmat, kernel_f, est_opt);
%
% Input
%    volmat: An array of volume matrices in the form N^3-by-N^3-by-K, where
%       N is the resolution and K is the number of matrices.
%    kernel_f: The Fourier transform of the convolution kernel in the form
%       M^3-by-M^3, where M >= N.
%    est_opt: The estimation options specifying the parameters of Ln. For
%       details, see the Sigma_est_opt argument to calc_Sigma_conv_cg.
%
% Output
%    volmat: The operator Ln applied to the input volmat, with added
%       regularization if desired.

function volmat = apply_conv_Ln(volmat, kernel_f, est_opt)
    est_opt = fill_struct(est_opt, ...
        'project', 0, ...
        'lambda', 0, ...
        'basis', [], ...
        'reweight_matrix', []);

    % TODO: Ensure basis is orthonormal.

    volmat = apply_conv(volmat, kernel_f, est_opt) + ...
        apply_reg(volmat, est_opt);
end

function volmat = apply_conv(volmat, kernel_f, est_opt)
    if ~isempty(est_opt.reweight_matrix)
        volmat = mat_conj(volmat, est_opt.reweight_matrix);
    end

    if ~isempty(est_opt.basis)
        volmat = mat_conj(volmat, est_opt.basis);
    end

    if est_opt.project
        N = size(vecmat_to_volmat(volmat), 1);

        basis_opt = legacy_basis(N-2, N, N);
        b2s3d = basis_to_space_3d(basis_opt);
    end

    if est_opt.project
        basis_project_Sigma_3d(volmat, b2s3d);
    end

    volmat = volmat_to_vecmat( ...
        conv_volmat(vecmat_to_volmat(volmat), ...
        vecmat_to_volmat(kernel_f)));

    if est_opt.project
        volmat = basis_project_Sigma_3d(volmat, b2s3d);
    end

    if ~isempty(est_opt.basis)
        volmat = mat_conj(volmat, est_opt.basis');
    end

    if ~isempty(est_opt.reweight_matrix)
        volmat = mat_conj(volmat, est_opt.reweight_matrix');
    end
end

function volmat = apply_reg(volmat, est_opt)
    if isnumeric(est_opt.lambda)
        volmat = est_opt.lambda*volmat;
    elseif isa(est_opt.lambda, 'function_handle')
        if ~isempty(est_opt.basis)
            volmat = mat_conj(volmat, est_opt.basis);
        end

        volmat = volmat_to_vecmat( ...
            volmat_apply_radial_filter( ...
            vecmat_to_volmat(volmat), est_opt.lambda));

        if ~isempty(est_opt.basis)
            volmat = mat_conj(volmat, est_opt.basis');
        end
    end
end

