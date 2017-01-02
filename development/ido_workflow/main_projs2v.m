
function v = main_projs2v (...
    SOURCE_FILE, BASE_PATH, N,...
    MOLEC_RADIUS, N_THETA, MAX_SHIFT, SHIFT_STEP,...
    VOTING_TICS_WIDTH, J_EIGS, J_WEIGHTS, REF_VOL,...
    S_WEIGHTS)
% CryoEM: reconstruct a 3D volume by its 2D projections images
% 
% Main wrapper written by Ido Greenberg, 2016


%% Setup

% Load Data
projs_uncleared = ReadMRC(SOURCE_FILE, 1, N);
n = size(projs_uncleared, 1);
log_message('Stage Done: Load Data');

% Configuration
% Necessary input: SOURCE_FILE, BASE_PATH, N
if ~exist('MOLEC_RADIUS','var') || isempty(MOLEC_RADIUS); MOLEC_RADIUS=0.4; end
if ~exist('N_THETA','var') || isempty(N_THETA); N_THETA=360; end
if ~exist('MAX_SHIFT','var') || isempty(MAX_SHIFT); MAX_SHIFT=floor(10*n/89); end
if ~exist('SHIFT_STEP','var') || isempty(SHIFT_STEP); SHIFT_STEP=1; end
if ~exist('VOTING_TICS_WIDTH','var') || isempty(VOTING_TICS_WIDTH); VOTING_TICS_WIDTH=1; end % GGG reconsider...
if ~exist('J_EIGS','var') || isempty(J_EIGS); J_EIGS=4; end
if ~exist('J_WEIGHTS','var') || isempty(J_WEIGHTS); J_WEIGHTS=true; end
if ~exist('REF_VOL','var') || isempty(REF_VOL); REF_VOL=''; end
if ~exist('S_WEIGHTS','var') || isempty(S_WEIGHTS); S_WEIGHTS=true; end

save(sprintf('%s/configuration', BASE_PATH),...
    'SOURCE_FILE', 'BASE_PATH', 'N', 'n',...
    'MOLEC_RADIUS', 'N_THETA', 'MAX_SHIFT', 'SHIFT_STEP',...
    'VOTING_TICS_WIDTH', 'J_EIGS', 'J_WEIGHTS', 'REF_VOL',...
    'S_WEIGHTS');
open_log(sprintf('%s/stdout.txt', BASE_PATH));
log_message('Source data: %s', SOURCE_FILE);
log_message('Output path: %s', BASE_PATH);
log_message('N=%d projection-images of size %dx%d', N,n,n);

% Clear Images
% Currently just applly circular mask
projs = mask_fuzzy(projs_uncleared, size(projs_uncleared,1)*MOLEC_RADIUS);
save(sprintf('%s/projections', BASE_PATH),...
    'projs_uncleared', 'projs');


%% Common Lines

% Polar Fourier Transform of projections
[npf, ~] = cryo_pft(projs, n, N_THETA, 'single');

% Common lines
[clstack,corrstack,~] =...
    cryo_clmatrix_gpu(npf, size(npf,3), 0, MAX_SHIFT, SHIFT_STEP, 0);
%     cryo_clmatrix_multi_gpu_opt(npf, size(npf,3), shifts.MAX_SHIFT, shifts.SHIFT_STEP, 5, n_cls, 1);

save(sprintf('%s/common_lines', BASE_PATH),...
    'clstack');
log_message('Stage Done: Common Lines');


%% Estimate Rotations

% Estimate relative rotations
[Rij0, r_valid_ks, r_good_ks, peak_width] = cryo_sync3n_estimate_all_Rijs...
    (clstack, N_THETA, VOTING_TICS_WIDTH);

save(sprintf('%s/relative_rotations', BASE_PATH),...
'Rij0');
log_message('Stage Done: Relative Rotations');

% J-synchronization
verbose = 2;
[J_sync,J_significance,J_eigs,J_sync_iterations,~] = ...
cryo_sync3n_Jsync_power_method (Rij0, J_EIGS, J_WEIGHTS, verbose);
Rij = cryo_sync3n_flip_handedness(J_sync, Rij0);

% Build 3NX3N Matrix S
S = cryo_sync3n_syncmatrix(Rij);

save(sprintf('%s/J_sync', BASE_PATH),...
    'J_sync','J_significance','J_eigs','J_sync_iterations','Rij','S');
log_message('Stage Done: J Synchronization');

% S Weighting
W = ones(N); % Default weights are all one (no weights)
if S_WEIGHTS
    % Estimate weights for the 3x3 blocks of S
    [W, Pij, scores_hist] = cryo_sync3n_syncmatrix_weights(Rij0);
end

save(sprintf('%s/weights', BASE_PATH),...
    'W','Pij','scores_hist');
log_message('Stage Done: Probabilistic Weighting');

% Estimate rotations from S
[rotations, S_eigs, ~] = cryo_sync3n_S_to_rot (S,10,W);

save(sprintf('%s/absolute_rotations', BASE_PATH),...
    'rotations','S_eigs');
log_message('Stage Done: Absolute Rotations');


%% Reconstruction

% Shifts
[shifts,~] = cryo_estimate_shifts(cryo_pft(projs,n,N_THETA,'single'),...
    rotations, MAX_SHIFT, SHIFT_STEP);

% Global Phase Flip
[projs, flipped] = cryo_globalphaseflip(projs);
if flipped; log_message('Note: all projections contrasts were flipped.\n'); end

save(sprintf('%s/shifts_and_phaseflip', BASE_PATH),...
    'shifts','flipped');
log_message('Stage Done: Shifts Estimation & Global Phase Flip');

% Molecule Reconstruction
[v, ~, ~, ~, ~, ~] = recon3d_firm( projs, rotations, -shifts, 1e-6, 30, zeros(n,n,n) );
if norm(imag(v(:)))/norm(v(:)) > 1e-3
    warning('Imaginary part of the reconstruction is not neglible: %f',...
        norm(imag(v(:)))/norm(v(:)));
end
v = real(v);

WriteMRC(v,1,sprintf('%s/vol.mrc', BASE_PATH));
log_message('Stage Done: Density Map Reconstruction');


%% Post Analysis

% Volume Alignment and FSC
if exist('REF_VOL','var') && ~isempty(REF_VOL)
    % Alignment
    v_ref = ReadMRC(REF_VOL);
    [~,~,v_aligned,~] = cryo_align_densities(v_ref,v);
    WriteMRC(v_aligned, 1, sprintf('%s/vol_ref-aligned.mrc',BASE_PATH));
    % FSC
    [res,fc] = fsc(v_ref, v_aligned);
    [~, fighandle] = plot_fsc(fsc);
    
    save(sprintf('%s/resolution', BASE_PATH),...
        'res','fc');
    log_message('Stage Done: Alignment vs. Reference');
    log_message('Resolution vs. Reference:\t%.2f [A]', res);
end

% Finish
memory_usage(whos);
log_message('Reconstruction Done.');

end
