function rotations=cryo_sync3n_estimate_rotations(clstack,L)
% CRYO_SYNC3N_ESTIMATE_ROTATIONS Estimate orientations from common lines.
%
% cryo_sync3n_estimate_rotations(clstack,L)
%   Estimate orientaitons corresponding to the common lines matrix given by
%   clstack. The angular rsolution used is L.
%   The function returns an array of dimensions 3x3xN, where N is the
%   number of images, where rotations(:,:,i)is the rotation corresponding
%   to the i'th image.
%
%   The rotations are estiamted using the 3Nx3N unweighted algorithm.
%
% Yoel Shkolnisky, March 2016.

% Paramters XXX scores_as_entries, verbose,

K=size(clstack,1);

log_message('Estimating pairwise rotations');
% Estimate the relative rotation Rij between each pair of images.
[Rijs, ~, ~, ~] = cryo_sync3n_estimate_all_Rijs(clstack, L);
        
n_eigs=3; % Number of eigenvalues of the J-sync matrix to compute.
scores_as_entries = 0; % Construct the J-sync matrix using only +1 and -1.
verbose=1;

log_message('Synchronizing relative orientations');
log_message('Computing %d eigenvalues of the J-synchronization matrix',n_eigs);
% Synchronize all relative rotations to a consistent handedness.
[J_sync,~,eigenvalues,itr,~] =...
    cryo_sync3n_Jsync_power_method(Rijs,n_eigs,scores_as_entries,verbose );
log_message('Converged on iterations %d',itr);
log_message('Eigenvalues of the J-synchronization matrix [ %s ]',num2str(eigenvalues.'));


% Flip the handedness of the estimates Rij according to hand estimated in
% J_sync. Estimates Rij from which J_sync==-1 are J-conjugated.
Rijs = cryo_sync3n_flip_handedness(J_sync, Rijs);

% Build 3KX3K synchronization matrix.
log_message('Building synchronization matrix');
S = cryo_sync3n_syncmatrix(Rijs);

% Estimate rotations from S.
[rotations, ~, ~] = cryo_sync3n_S_to_rot (S, 10, -1);