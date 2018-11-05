% This script illustrates the basic ab initio reconstruction functionality of
% the ASPIRE toolbox on simulated data, using common lines to recover the
% viewing angles from clean images and reconstructing the volume using the
% least-squares volume estimator.

%%% Parameters %%%

L = 33;                 % Size of images.
n = 128;                % Number of images.

SNR = 32;               % Signal-to-noise ratio of images.

n_r = ceil(L/2);        % Number of radial nodes in polar Fourier transform.
n_theta = 36;           % Number of angular nodes in polar Fourier transform.

%%% Simulate data %%%

% Load the 'cleanrib' volume, corresponding to the experimentally obtained EM
% map of a 50S ribosome.
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

%%% Estimate viewing angles %%%

% Calculate polar Fourier transforms of images.
pf = cryo_pft(ims, n_r, n_theta);

% Estimate common-lines matrix.
clstack_est = cryo_clmatrix(pf);

% Calculate the "true" common-lines matrix from true rotations.
clstack_true = clmatrix_cheat(rots_true, n_theta);

% Construct syncronization matrix from common lines.
S_est = cryo_syncmatrix_vote(clstack_est, n_theta);

% Estimate rotations using synchronization.
inv_rots_est = cryo_syncrotations(S_est);

% Align our estimates rotations to the true ones. This lets us compute the MSE,
% but also allows us to more easily compare the reconstructed volume.
inv_rots_aligned = align_rots(inv_rots_est, inv_rots_true);

%%% Estimate volume %%%

% Set up parameters for volume estimation.
params = struct();
params.rot_matrices = inv_rots_aligned;     % Estimated rotations.
params.ctf = ones(L*ones(1, 2));            % CTFs (none here).
params.ctf_idx = ones(1, n);                % CTF indices (all the same).
params.ampl = ones(1, n);                   % Amplitude multipliers (all one).
params.shifts = zeros(2, n);                % Shifts (none here).

% Set up basis in which to estimate volume. Here, it is just the standard Dirac
% basis, where each voxel is a basis vector.
basis = ffb_basis(L*ones(1, 3));

% Set up options for the volume estimation algorithm.
mean_est_opt = struct();
mean_est_opt.verbose = false;               % Don't output progress info.
mean_est_opt.max_iter = 10;                 % Maximum number of iterations.
mean_est_opt.rel_tolerance = 1e-3;          % Stopping tolerance.
mean_est_opt.half_pixel = true;             % Center volumes around half pixel.
mean_est_opt.verbose = 1;                   % Print progress information.

% Estimate volume using least squares.
vol_est = cryo_estimate_mean(ims, params, basis, mean_est_opt);

%%% Evaluate results %%%

% Calculate proportion of correctly identified common lines.
cl_prop = comparecl(clstack_est, clstack_true, n_theta, 10);

% Calculate MSE of rotations.
mse_rots = tnorm(inv_rots_aligned-inv_rots_true)^2/n;

% Calculate relative MSE of volume estimation.
nrmse_vol = tnorm(vol_est-vol_true)/tnorm(vol_true);

% Reproject estimated volume and compare with clean images.
ims_est = cryo_project(vol_est, rots_true);
ims_est = permute(ims_est, [2 1 3]);

nrmse_ims = tnorm(ims_est-ims_clean)/tnorm(ims_clean);

fprintf('%-40s%20g\n', 'Proportion of correct common lines:', cl_prop);
fprintf('%-40s%20g\n', 'MSE of rotations:', mse_rots);
fprintf('%-40s%20g\n', 'Estimated volume normalized RMSE:', nrmse_vol);
fprintf('%-40s%20g\n', 'Reprojected images normalized RMSE:', nrmse_ims);
