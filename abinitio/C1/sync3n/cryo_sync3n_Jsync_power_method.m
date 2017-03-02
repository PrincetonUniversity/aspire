function [J_sync,J_significance,eigenvalues,itr,dd] =...
    cryo_sync3n_Jsync_power_method (Rij , n_eigs , scores_as_entries , verbose , measure_time )
%SIGNS_POWER_METHOD J-synchronization.

%cryo_sync3n_Jsync_power_method (Rij , n_eigs , scores_as_entries , verbose , measure_time )
%   Calculate n_eigs eigenvalues & eigenvectors of the J-synchronization
%   matrix using power method.
%
% As the J-synchronization matrix is of size (N choose 2)x(N choose 2), the
% matrix uses the power method of compute the eigenvalues and eigenvectors,
% while constructing the matrix on-the-fly.
% The multiplication of the J-synchronization matrix by a vector is
% implemented in signs_times_v_mex().
%
% Input:
% Rij   3X3X(N choose 2) array of estimates of Ri'Rj, computed, for
%       example, by cryo_sync3n_estimate_all_Rijs.
% n_eigs    Number of eigenvalues and eigenvectors to compute.
% score_as_entries  Whether to use triplets scores as the matrix
%       entries or just the signs.
% verbose   0 for silent, 1 for summary printouts, 2 for detailed
%       printouts.
%
% Output:
% J_sync    (N choose 2) dimensional vector with the k'th entry indicating
%       whether Rij(:,:,k) should be J-conjugated or not to acheive global
%       handedness consistency. The vector consists only of +1 and -1.
% J_significance    (N choose 2) dimensional vector with the k'th entry
%       indicating the reliability of the corresponding entry in J_sync,
%       computed based on the eigenvector of the signs matrix. This is in
%       fact the absolute value of the top eigevector of the
%       J-synchronization matrix.
% eigenvalues   Top n_eigs eigenvalues of the J_synchronization matrix.
% itr   Number of iterations required for the power method to converege.
% GGG Yoel's question: How did you determine convergence?
% dd    Magnitude of last change in the estimated eigenvector.
%
%
% Written by Ido Greenberg, 2015.
% Revised by Yoel Shkolnisky, March 2016.

N=(1+sqrt(1+8*size(Rij,3)))/2; % Extract N from the number of relative rotaitons in Rij.
% The number of relative rotation for a given N is (N choose 2) so,
% solve for N from this number.
assert(N==round(N));  % Make sure we got an integer.


% input validation
if (n_eigs <= 0); error('n_eigs must be positive integer!'); end
if ~exist('verbose','var'); verbose=0; end
if ~exist('measure_time','var'); measure_time=false; end
pairs_scores = []; % currently disabled. intended to allow certain weighting of the J-sync matrix.

% constants
epsilon = 1e-2;
MAX_ITERATIONS = 100;
if verbose >= 2
    log_message('Power method for signs matrix started, up to eigenvector\naccuracy goal of %f, and limited by %d iterations.',...
        epsilon, MAX_ITERATIONS);
end

% set input to fit the mex interface
if scores_as_entries
    scores_as_entries = 1;
else
    scores_as_entries = 0;
end

% initialization
N_pairs = size(Rij,3);
vec = rand(N_pairs,n_eigs);
[vec,~] = qr(vec,0);
dd = 1;
itr = 0;

% power method iterations
if measure_time; t = toc; end
while itr < MAX_ITERATIONS && dd > epsilon
    itr = itr + 1;
    [vec_new,n_cores] = signs_times_v_mex(N,Rij,vec,n_eigs,scores_as_entries,pairs_scores,0);
    vec_new = reshape(vec_new,size(vec));
    [vec_new,eigenvalues] = qr(vec_new,0);
    dd = norm(vec_new(:,1)-vec(:,1));
    vec = vec_new;
    if verbose >= 2
        log_message('Iteration %02d (%02d cores used):\t||curr_evec-last_evec|| = %.3f',...
            itr, n_cores, dd);
    end
end
if measure_time; t=(toc-t)/60; end

% set output format
eigenvalues = diag(eigenvalues);
J_significance = abs(vec(:,1));
J_sync = sign(vec(:,1)); % only signs of the eigenvector entries are needed


% print eigenvalues & number of iterations needed for the power method
if verbose >= 1
    log_message('Outer J Synchronization:\n\titerations needed: %d', itr);
    if measure_time; log_message('\ttime consumed: %f [min]', t); end
    log_message('\tsigns matrix eigenvalues (in practice vs. theoretical):');
    log_message('Note: theoretical values are in the noiseless case with no weights (scores_as_entries=0).');
    log_message('\t\t%.0f\t(%.0f)', eigenvalues(1), 2*N-4);
    for i=2:min(n_eigs,5); log_message('\t\t%.0f\t(%.0f)', eigenvalues(i), N-4); end
end

end
