% This script illustrates the basic reconstruction functionality of the ASPIRE
% toolbox on simulated data using the true viewing angles.

%%% Parameters %%%

L = 33;                 % Size of images.
n = 128;                % Number of images.

SNR = 32;               % Signal-to-noise ratio of images.

%%% Simulate data %%%

% Load the 'cleanrib' volume, corresponding to the experimentally obtained EM
% map of a 70S ribosome.
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, L*ones(1, 3));

% Generate random rotations. These are the true viewing directions.
rots_true = rand_rots(n);
inv_rots_true = invert_rots(rots_true);

% Generate clean image by projecting volume along rotations.
ims_clean = cryo_project(vol_true, rots_true);

% `cryo_project` does not generate images compatible with the other functions in
% the package, so we need to switch the x and y axes.
ims_clean = permute(ims_clean, [2 1 3]);

% Add noise at desired SNR.
power_clean = tnorm(ims_clean).^2/numel(ims_clean);
sigma_noise = sqrt(power_clean/SNR);
ims = ims_clean + sigma_noise*randn(size(ims_clean));

%%% Estimate volume %%%

% Set up parameters for volume estimation.
params = struct();
params.rot_matrices = inv_rots_true;        % True viewing angles.
params.ctf = ones(L*ones(1, 2));            % CTFs (none here).
params.ctf_idx = ones(1, n);                % CTF indices (all the same).
params.ampl = ones(1, n);                   % Amplitude multipliers (all one).
params.shifts = zeros(2, n);                % Shifts (none here).

% Set up basis in which to estimate volume. Here, it is just the standard Dirac
% basis, where each voxel is a basis vector.
basis = dirac_basis(L*ones(1, 3));

% Set up options for the volume estimation algorithm.
mean_est_opt = struct();
mean_est_opt.verbose = false;               % Don't output progress info.
mean_est_opt.max_iter = 10;                 % Maximum number of iterations.
mean_est_opt.rel_tolerance = 1e-3;          % Stopping tolerance.
mean_est_opt.half_pixel = true;             % Center volumes around half pixel.

% Estimate volume using least squares.
vol_est = cryo_estimate_mean(ims, params, basis, mean_est_opt);

%%% Evaluate results %%%

% Calculate relative MSE of volume estimation.
nrmse_vol = tnorm(vol_est-vol_true)/tnorm(vol_true);

% Reproject estimated volume and compare with clean images.
ims_est = cryo_project(vol_est, rots_true);
ims_est = permute(ims_est, [2 1 3]);

nrmse_ims = tnorm(ims_est-ims_clean)/tnorm(ims_clean);

fprintf('%-40s%20g\n', 'Estimated volume normalized RMSE:', nrmse_vol);
fprintf('%-40s%20g\n', 'Reprojected images normalized MSE:', nrmse_ims);
