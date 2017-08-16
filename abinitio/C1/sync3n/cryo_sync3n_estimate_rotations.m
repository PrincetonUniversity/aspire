function [rotations,S]=cryo_sync3n_estimate_rotations(clstack,L,use_J_weights,use_S_weights)
% CRYO_SYNC3N_ESTIMATE_ROTATIONS Estimate orientations from common lines.
%
% cryo_sync3n_estimate_rotations(clstack,L)
%   Estimate orientaitons corresponding to the common lines matrix given by
%   clstack. The angular rsolution used is L.
%   The function returns an array of dimensions 3x3xN, where N is the
%   number of images, where rotations(:,:,i)is the rotation corresponding
%   to the i'th image. The rotations are estiamted using the 3Nx3N
%   unweighted algorithm.
%
% cryo_sync3n_estimate_rotations(clstack,L,use_weights)
%   Weight the J-synchronization matrix as well as the 3N synchronization
%   matrix to improve robustness. Default is use_weights=0 (no weights).
%
% Yoel Shkolnisky, April 2016.

if nargin<3
    use_J_weights=0;
end

if nargin<4
    use_S_weights=0;
end

K=size(clstack,1);

log_message('Estimating pairwise rotations');
% Estimate the relative rotation Rij between each pair of images.
[Rijs, ~, ~, ~] = cryo_sync3n_estimate_all_Rijs(clstack, L);
        
% If use_weights=0, then construct the J-sync matrix using only +1 and -1.
if use_J_weights
   use_J_weights=1;
end

n_eigs=3; % Number of eigenvalues of the J-sync matrix to compute.
verbose=1;

log_message('Synchronizing relative orientations');
log_message('Computing %d eigenvalues of the J-synchronization matrix',n_eigs);
% Synchronize all relative rotations to a consistent handedness.
[J_sync,scores,eigenvalues,itr,~] =...
    cryo_sync3n_Jsync_power_method(Rijs,n_eigs,use_J_weights,verbose );
log_message('Converged on iteration %d',itr);
log_message('Eigenvalues of the J-synchronization matrix [ %s ]',num2str(eigenvalues.'));

% Flip the handedness of the estimates Rij according to hand estimated in
% J_sync. Estimates Rij from which J_sync==-1 are J-conjugated.
Rijs = cryo_sync3n_flip_handedness(J_sync, Rijs);

% Build 3KX3K synchronization matrix.
log_message('Building synchronization matrix');
S = cryo_sync3n_syncmatrix(Rijs);

% Simulations show that the best-performing weights for the blocks of S are
% the values of the top eigenvector of the J-synchronization matrix.
% These values are stored in the variable scores returned by the function
% cryo_sync3n_Jsync_power_method above.
% Update: the probabilities of the relative rotations {Rij} to be indicative ("good")
% can be estimated by the significance of the local J-synchronization between every
% triplet Rij,Rjk,Rik. Those probabilities turn out to have better performance as weights,
% both analytically and empirically.
W=ones(K); % Default weights are all one (no weights)
if use_S_weights
    % Estimate weights for the 3x3 blocks of S
    W=cryo_sync3n_syncmatrix_weights(Rijs);
end

% Estimate rotations from S.
[rotations, ~, ~] = cryo_sync3n_S_to_rot (S, 10,W,-1);
