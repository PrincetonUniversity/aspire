
function v = main_projs2v (...
    SOURCE_PATH, BASE_PATH, S_WEIGHTS, REF_VOL)
% CryoEM: reconstruct a 3D volume by its 2D projections images
% 
% Main wrapper written by Ido Greenberg, 2016


%% Setup
% Load beginning of reconstruction from external source

% Load configuration
% Necessary input: SOURCE_PATH, BASE_PATH, S_WEIGHTS
conf = load(sprintf('%s/configuration', SOURCE_PATH));
N = conf.N;
MOLEC_RADIUS = conf.MOLEC_RADIUS;
N_THETA = conf.N_THETA;
MAX_SHIFT = conf.MAX_SHIFT;
SHIFT_STEP = conf.SHIFT_STEP;
VOTING_TICS_WIDTH = conf.VOTING_TICS_WIDTH;
J_EIGS = conf.J_EIGS;
J_WEIGHTS = conf.J_WEIGHTS;
if ~exist('REF_VOL','var') || isempty(REF_VOL); REF_VOL=conf.REF_VOL; end
n = conf.n;

save(sprintf('%s/configuration', BASE_PATH),...
    'SOURCE_PATH', 'BASE_PATH', 'N', 'n',...
    'MOLEC_RADIUS', 'N_THETA', 'MAX_SHIFT', 'SHIFT_STEP',...
    'VOTING_TICS_WIDTH', 'J_EIGS', 'J_WEIGHTS', 'REF_VOL',...
    'S_WEIGHTS');

% Load the relevant results from the reconstruction until Rij Synchronization
projs = load(sprintf('%s/projections', SOURCE_PATH));
projs = projs.projs;
Rij0 = load(sprintf('%s/relative_rotations', SOURCE_PATH));
Rij0 = Rij0.Rij0;
S = load(sprintf('%s/J_sync', SOURCE_PATH));
S = S.S;

save(sprintf('%s/restored_from_source', BASE_PATH),...
    'projs', 'Rij0', 'S');

% Print initial logs
open_log(sprintf('%s/stdout.txt', BASE_PATH));
log_message('\nBeginning reconstruction...');
log_message('Source for Synced {Rij}: %s', SOURCE_PATH);
log_message('Output path: %s', BASE_PATH);
log_message('N=%d projection-images of size %dx%d', N,n,n);
log_flush();


%% Estimate Rotations

% S Weighting
W = ones(N); % Default weights are all one (no weights)
if S_WEIGHTS
    % Estimate weights for the 3x3 blocks of S
    [W, Pij, scores_hist] = cryo_sync3n_syncmatrix_weights(Rij0);
end

save(sprintf('%s/weights', BASE_PATH),...
    'W','Pij','scores_hist');
log_message('Stage Done: Probabilistic Weighting');
log_flush();

% Estimate rotations from S
[rotations, S_eigs, ~] = cryo_sync3n_S_to_rot (S,10,W);

save(sprintf('%s/absolute_rotations', BASE_PATH),...
    'rotations','S_eigs');
log_message('Stage Done: Absolute Rotations');
log_flush();


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
log_flush();

% Molecule Reconstruction
[v, ~, ~, ~, ~, ~] = recon3d_firm( projs, rotations, -shifts, 1e-6, 30, zeros(n,n,n) );
if norm(imag(v(:)))/norm(v(:)) > 1e-3
    warning('Imaginary part of the reconstruction is not neglible: %f',...
        norm(imag(v(:)))/norm(v(:)));
end
v = real(v);

WriteMRC(v,1,sprintf('%s/vol.mrc', BASE_PATH));
log_message('Stage Done: Density Map Reconstruction');
log_flush();


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
close_log();

end
