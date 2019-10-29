% This script illustrates the covariance Wiener filtering functionality of the
% ASPIRE toolbox, implemented by estimating the covariance of the unfiltered
% images in a Fourier-Bessel basis and applying the Wiener filter induced by
% that covariance matrix.

%%% Parameters %%%

L = 8;                 % Size of images
n = 32;               % Number of images n
% n = 1024;               % Number of images n
SNR = 1;                % Signal-to-noise ratio of images.

pixel_size = 5;         % Pixel size of the images (in angstroms).

defocus_min = 1.5e4;    % Minimum defocus value (in angstroms).
defocus_max = 2.5e4;    % Maximum defocus value (in angstroms).
defocus_ct = 7;         % Number of defocus groups.

basis = ffb_basis(L*ones(1, 2), L, 1);

% Initialize the random number generator.
initstate(0);

%%% Simulation setup %%%

% Set up the arrays to hold the CTFs, both as 2D filters and as block diagonal
% matrices in the Fourier-Bessel basis.
h = zeros([L*ones(1, 2) defocus_ct]);
h_fb = cell(1, defocus_ct);

% Set up the default CTF parameters. The only one we vary is the defocus.
ctf_params = struct();
ctf_params.voltage = 200;
ctf_params.spherical_aberration = 2;
ctf_params.amplitude_contrast = 0.1;

% Calculate the defocus values which will define the CTFs.
defocus = linspace(defocus_min, defocus_max, defocus_ct);

% For each defocus, generate a function handle representing the CTF as a
% function of radial frequency, then evaluate that on a grid and calculate
% its block diagonal representation in the Fourier-Bessel basis.
for k = 1:numel(defocus)
    ctf_params.defocus = defocus(k);

    h_fun = cryo_radial_ctf_handle(pixel_size, ctf_params);

    h(:,:,k) = radial_fun_to_image(h_fun, L);
    h_fb{k} = radial_filter_to_fb_mat(h_fun, basis);
end

% Randomly assign a CTF to each image.
h_idx = randi(numel(defocus), [1 n]);

% Load the 'cleanrib' volume, corresponding to the experimentally obtained EM
% map of a 70S ribosome.
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol = cryo_downsample(f.volref, L*ones(1, 3));

% Generate random rotations and project.
initstate(0);
rots = rand_rots(n);
proj_clean = cryo_project(vol, rots);

% Prepare an array to hold the filtered projections.
proj_ctf_clean = zeros(size(proj_clean), class(proj_clean));

% For each defocus group, find the images that are assigned to that CTF, and
% filter them.
for k = 1:numel(defocus)
    mask = find(h_idx == k);

    proj_ctf_clean(:,:,mask) = im_filter(proj_clean(:,:,mask), h(:,:,k));
end

% Add noise at the desired SNR to the filtered images to get our final
% projections.
power_clean = tnorm(proj_ctf_clean)^2/numel(proj_ctf_clean);
noise_var = power_clean/SNR;
initstate(0);
proj = proj_ctf_clean + sqrt(noise_var)*randn2(size(proj_ctf_clean));

%%% Estimation %%%

% Expand the images, both clean and noisy, in the Fourier-Bessel basis. This
% can be done exactly (that is, up to numerical precision) using the
% `basis.expand` function, but for our purposes, an approximation will do.
% Since the basis is close to orthonormal, we may approximate the exact
% expansion by applying the adjoint of the evaluation mapping using 
% `basis.evaluate_t`.
coeff_clean = basis.evaluate_t(proj_clean);
coeff = basis.evaluate_t(proj);

% Given the clean Fourier-Bessel coefficients, we can estimate the symmetric
% mean and covariance. Note that these are not the same as the sample mean and
% covariance, since these functions use the rotational and reflectional
% symmetries of the distribution to constrain to further constrain the
% estimate. Note that the covariance matrix estimate is not a full matrix,
% but is block diagonal. This form is a consequence of the symmetry
% constraints, so to reduce space, only the diagonal blocks are stored. The
% mean and covariance estimates will allow us to evaluate the mean and
% covariance estimates from the filtered, noisy data, later.
mean_coeff = fb_rot_mean(coeff_clean, basis);
covar_coeff = fb_rot_covar(coeff_clean, mean_coeff, basis);

% We now estimate the mean and covariance from the Fourier-Bessel
% coefficients of the noisy, filtered images. These functions take into
% account the filters applied to each image to undo their effect on the
% estimates. For the covariance estimation, the additional information of
% the estimated mean and the variance of the noise are needed. Again, the
% covariance matrix estimate is provided in block diagonal form.
mean_coeff_est = fb_rot_mean_ctf(coeff, h_fb, h_idx, basis);
% covar_coeff_est = fb_rot_covar_ctf(coeff, h_fb, h_idx, ...
%    mean_coeff_est, noise_var, basis);

covar_est_opt = struct();
covar_est_opt = fill_struct(covar_est_opt, ...
        'shrinker', 'frobenius_norm', ...
        'verbose', 0, ...
        'max_iter', 250, ...
        'rel_tolerance', 1e-12);
covar_coeff_est = fb_rot_covar_ctf(coeff, h_fb, h_idx, ...
    mean_coeff_est, noise_var, basis, covar_est_opt);

% Estimate the Fourier-Bessel coefficients of the underlying images using a
% Wiener filter. This Wiener filter is calculated from the estimated mean,
% covariance, and the variance of the noise. The resulting estimator has
% the lowest expected mean square error out of all linear estimators.
coeff_est = fb_rot_wiener_ctf(coeff, h_fb, h_idx, mean_coeff_est, ...
    covar_coeff_est, noise_var, basis);

% Convert Fourier-Bessel coefficients back into 2D images.
proj_est = basis.evaluate(coeff_est);

%%% Evaluate results %%%

% Calculate the difference between the estimated covariance and the "true"
% covariance estimated from the clean Fourier-Bessel coefficients.
covar_coeff_diff = blk_diag_add(covar_coeff, ...
    blk_diag_mult(-1, covar_coeff_est));

% Calculate the deviation between the clean estimates and those obtained from
% the noisy, filtered images.
diff_mean = norm(mean_coeff_est-mean_coeff)/norm(mean_coeff);
diff_covar = blk_diag_norm(covar_coeff_diff)/blk_diag_norm(covar_coeff);

% Calculate the normalized RMSE of the estimated images.
nrmse_ims = tnorm(proj_est-proj_clean)/tnorm(proj_clean);

fprintf('%-50s%20g\n', 'Deviation of the noisy mean estimate:', diff_mean);
fprintf('%-50s%20g\n', 'Deviation of the noisy covariance estimate:', diff_covar);
fprintf('%-50s%20g\n', 'Estimated images normalized RMSE:', nrmse_ims);
