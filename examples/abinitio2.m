% Basic abinitio reconstruction - example 2.
%
% The script demonstrates calling various functions of the aspire package.
%
% Follow comments below.
%
% Yoel Shkolnisky, September 2013.

clear;

%% Generate simulated projections
% Generate 200 simulated projections of size 65x65.
% For simplicity, the projections are centered.
n=65;
K=200;
SNR=1/4;
log_message('Generating %d simulated noisy projections of size %dx%d',K,n,n);
%SNR=1000; % No noise
[projs,noisy_projs,~,rots_ref]=cryo_gen_projections(n,K,SNR);
viewstack(noisy_projs,5,5) % Show some noisy projections

% Find reference common lines 
n_theta=360;
log_message('Computing ground-truth common lines');
[ref_clstack,~]=clmatrix_cheat(rots_ref,n_theta);

% Estimate common lines from noisy projections
log_message('Estimating common lines between noisy projections');
masked_projs=mask_fuzzy(noisy_projs,33); % Applly circular mask
viewstack(masked_projs,5,5) % Show some noisy projections

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_r=65;
[npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

% Find common lines from projections
max_shift=0;
shift_step=1;
clstack = cryo_clmatrix(npf,-1,1,max_shift,shift_step);
prop=comparecl( clstack, ref_clstack, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);

%% Estimate orientations using sychronization.
log_message('Estimating orientations of projections using synchronization'); 
S=cryo_syncmatrix_vote(clstack,n_theta);
[rotations,diff,mse]=cryo_syncrotations(S,rots_ref);
fprintf('MSE of the estimated rotations: %f\n\n',mse);
%fprintf('MSE of the estimated rotations: %f\n\n', ...
%   check_MSE(rotations,rots_ref));


%% 3D inversion
log_message('Reconstructing 3D density from projections');
params = struct();
params.rot_matrices = rotations;
params.ctf = ones(size(noisy_projs, 1)*ones(1, 2));
params.ctf_idx = ones(size(noisy_projs, 3), 1);
params.shifts = zeros(2, size(noisy_projs, 3));
params.ampl = ones(size(noisy_projs, 3), 1);

basis = dirac_basis(size(noisy_projs, 1)*ones(1, 3));

v = cryo_estimate_mean(noisy_projs, params, basis);

fname='example2.mrc';
WriteMRC(v,1,fname); % Output density map reconstructed from projections.
log_message('Reconstructed density saved to %s',fname);
log_message('Done!');