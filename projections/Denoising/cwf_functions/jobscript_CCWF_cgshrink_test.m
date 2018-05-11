% JOBSCRIPT_CCWF_CGSHRINK_TEST Perform CWF denoising
%
% Usage
%    [denoised_coeff, skip_flags] = jobscript_CCWF_cgshrink_test(index, ...
%       CTF_rad_all, basis, sample_points, mean_coeff, coeff_pos_k, noise_v);
%
% Input
%    index: A vector of n indices corresponding to the CTF group of each image.
%    CTF_rad_all: The CTFs evaluated along the radial sample points.
%    basis: The Fourier-Bessel basis functions. Obtained from precomp_fb.
%    sample_points: The sample points for the basis. Obtained from precomp_fb.
%    mean_coeff: The coefficients of the mean image for k = 0.
%    coeff_pos_k: The Fourier-Bessel coefficients of the images. This is a
%       cell array of the same format as the output of FBcoeff_nfft.
%    noise_v: The noise variance of hte images.
%
% Output
%    denoised_coeff: A cell array of the denoised coefficients, with each cell
%       containing the coefficients for a different angular frequency, like
%       coeff_pos_k.
%    skip_flags: An array containing the skip flags for each angular frequency,
%       indicating for which frequencies the CG step was skipped.

% Written by Tejal Bhamre - Oct 2015
% Reformatted, refactored, and documented by Joakim Anden - 2018-Apr-19

function [denoised_coeff, skip_flags] = jobscript_CCWF_cgshrink_test(index, CTF_rad_all, basis, sample_points, mean_coeff, coeff_pos_k, noise_v)
    nim = length(index);
    ang_freqs = basis.ang_freqs;

    for i = 1:max(index)
        weight(i) = length(find(index==i));
    end

    denoised_coeff = cell(length(coeff_pos_k), 1);

    % Only pos k, to be compatibel with Jane's ffb
    for k = unique(ang_freqs)'
        if k> = 0
            tmp = coeff_pos_k{k+1};

            for i = 1:max(index)
                DD{i} = real(tmp(:,find(index==i))*tmp(:,find(index==i))'/length(find(index==i)));

                A{i} = calc_fb_CTF(CTF_rad_all(:,i), basis.Phi_ns{k+1}, sample_points);
            end

            % No regularization for white noise case
            regu = 0;

            [C, ~, ~, ~, skip_flag] = solve_cg_shrink(A, DD, weight, noise_v, nim, k, regu);

            skip_flags(k+1) = skip_flag;

            % Ensure positivity of the estimated covariance matrix
            [U, S] = eig(C);
            C = U*diag(max(diag(S), 0))*U';

            for i = 1:max(index)
                h_posk = C*A{i}'*inv(A{i}*C*A{i}'+noise_v*eye(size(A{i}, 1)));

                if (k == 0)
                    denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) + repmat(mean_coeff{k+1}, 1, numel(find(index==i)));
                else
                    denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i));
                end
            end
        end
    end
end
