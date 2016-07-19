function W = cryo_sync3n_syncmatrix_weights(scores)
%
% CRYO_SYNC3N_SYNCMATRIX_WEIGHTS   Determine weights for 3Nx3N synchronization matrix.
%
% cryo_sync3n_syncmatrix_weights(scores)
%   Determine weights matrix for the blocks of the 3N X 3N matrix S
%   (which are the relative rotations Rij between the projections
%   directions). The weights are used to improve the accuracy of the
%   estimated fotations, by giving more weight to reliable relative
%   rotations.
%
% Input:
% scores    A list of the weight to assign to each relative rotation in the
%           matrix S.  If we want to estimate N rotations, then scores is a
%           vector of length N choose 2.
%
% Output:
% W         NXN symmetric weights matrix without rows-normalization. 
%
% See cryo_sync3n_estimate_rotations for a useage example. 
%
% This function was revised from 
% set_weights...
%    set_weights(N, component_indices, W0, method, scores, stack, keep_initial_weights, iteration, N_iterations, w_normalization)
%
% Revised by Yoel Shkolnisky, April 2016, based on the function set_weights
% by Ido Greenberg, 2015 


N=(1+sqrt(1+8*numel(scores)))/2; % Extract N from the number of scores.
    % The number of scores for a given N is (N choose 2) as each pair of
    % rotations Rij (i<j) has a reliability score. Solve for N from this
    % number.
assert(N==round(N));  % Make sure we got an integer.

W=eye(N);
scores=scores(:);
scores = scores./max(scores);
idx = 0; % pair index
for i = 1:N
    for j = (i+1):N
        idx = idx + 1;
        W(i,j) = scores(idx);
        W(j,i) = scores(idx);
    end
end

end

