function [R, evalues, orthogonality_error] = cryo_sync3n_S_to_rot... 
    (S, n_eigs, W, ORTHOGONALITY_THRESHOLD)
% CRYO_SYNC3N_S_TO_ROT  Extract rotations from synchronization matrix.
%
% cryo_sync3n_S_to_rot (S, n_eigs, W, ORTHOGONALITY_THRESHOLD)
%   Compute rotations from the synchronization matrix S. Given a 3Nx3N
%   synchronization matrix S, whose (i,j)'th block of size 3x3 is an
%   estimate Ri'*Rj, for estimate the rortaions Ri. The function
%   computes n_eigs of the matrix S. n_eigs must be at least 3. If n_eigs
%   is larger, the fourth eigenvalue and on are used only to display the
%   gap in the spectrum of S, and are not used for computing the rotations.
%   The blocks of the matrix S are weighted according to the weights in W.
%   If S is of size 3Nx3N, then W is of size NxN whose (i,j) entry is the 
%   reliability of the (i,j)'th 3x3 block of S. 
%   When S is consturcted from noisy data, its top three eigenvectors are
%   not exactly the requested rotations. ORTHOGONALITY_THRESHOLD is used to
%   determine which estimated matrices are not orthogonal. Each estimated
%   matrix is considered orthogonal if the difference from its closest
%   orthogonal matrix (in Frobenius norm) is less than
%   ORTHOGONALITY_THRESHOLD. The function reports how many estimated
%   matrices are not orthogonal. Set ORTHOGONALITY_THRESHOLD to a negative
%   value to bypass this test. In any case, all returned matrices are
%   orthogonalized by projecting each estimated matrix to its nearest
%   orthogonal matrix.
%
% cryo_sync3n_S_to_rot (S, n_eigs, W)
%   Use no ORTHOGONALITY_THRESHOLD.
%
% cryo_sync3n_S_to_rot (S, n_eigs, ORTHOGONALITY_THRESHOLD)
%   No weights are applied to the blocks of the matrix S.
%
% cryo_sync3n_S_to_rot (S, n_eigs)
%   Use no ORTHOGONALITY_THRESHOLD and not weights.
%
% cryo_sync3n_S_to_rot (S)
%   Use neigs=3, no ORTHOGONALITY_THRESHOLD and not weights.
%
% Input:
%   S   3NX3N Synchronization matrix whose blocks are the relative
%       rotations, as returned from cryo_sync3n_syncmatrix.
%   n_eigs  Number of eigenvalues to compute
%   W   NXN weights matrix (without rows normalization).
%   ORTHOGONALITY_THRESHOLD    The function estimates the rotations from
%       the top three eigenvectors of the matrix S. In the noiseless case,
%       theses eigenvectors are exactly the rotations. In the presence of
%       noise, the resulting 3x3 matrices are not orthogonal, but we make
%       them orthogonal before returning them in R. The function counts the
%       number of matrices that deviate from their orthogonalized form (in
%       Frobenius norm) by more than ORTHOGONALITY_THRESHOLD.
%
% Output:
%   R   3X3XN rotations matrices.
%   eigenvalues The top n_eigs eigenvalues of the normalized, weighted S.
%   orthogonality_error Array of length N of the distances
%                       ||Ri-orthogonalize(Ri)||_F. 
%
% This function was revised from
%   rotations (S, W, Dhalf, N, n_eigs, ORTHOGONALITY_THRESHOLD)
% Parameter N was removed and is calculated from the dimensions of S.
%
% Yoel Shkolnisky, March 2016.
% Based on code by Ido Greenberg, 2015

%(S, n_eigs, W, ORTHOGONALITY_THRESHOLD)

N = size(S,1)/3;

if ~exist('n_eigs','var'); n_eigs=3; end
if ~exist('W','var'); W=ones(N); end
if ~exist('ORTHOGONALITY_THRESHOLD','var'); ORTHOGONALITY_THRESHOLD=1; end % default: ignore

% Input validation
if n_eigs < 3
    error('n_eigs must be >= 3 integer! otherwise H cannot be deduced from Sync Matrix. Currently n_eigs = %d', n_eigs);
end

% Prepare weights for the blocks of the matrix S.
D = mean(W,2); % We use mean and not sum so that after normalizing W the 
               % sum of its rows will be N, as required.
nulls = abs(D)<1.0e-13;
Dhalf = D;
Dhalf(~nulls) = Dhalf(~nulls).^-0.5;
Dhalf(nulls) = 0;
Dhalf = diag(Dhalf);

% Assert Rows Normalization
W_normalized = (Dhalf.^2)*W;
nzidx = sum(W_normalized,2) ~= 0; % do not consider erased rows
err = norm(sum(W_normalized(nzidx,:),2)-N); % GGG maybe max(sum) should be taken?
if err > 1e-10
    warning('Large Weights Matrix Normalization Error: %f', err);
end
%assert(norm(sum(W_normalized(nzidx,:),2)-N)<1.0e-9, 'Weights Matrix Normalization Error');

W = kron(W,ones(3));  % Make W of size 3Nx3N
Dhalf = diag(kron(diag(Dhalf),ones(3,1))); % Make Dhalf of size 3Nx3N

% Compute eigenvectors
[evecs, evalues] = eigs( Dhalf * (W.*S) * Dhalf , n_eigs);

% eigs() returns eigenvalues sorted by absolute value. we wish to sort
% the first 3 by real value, and the rest by either real or absolute value.
% the choice between both cases is used only for computing the spectral gap,
% and not for the reconstruction itself.
SORT_EIGS_BY_ABS_VALUE = false; % false = dont count large negative eigenvalues.
evalues = diag(evalues);
[~,ids] = sort(evalues, 'descend'); % (assert 3 largest eigenvalues are first)
evalues = evalues(ids);
evecs = evecs(:,ids);
if n_eigs > 3
    if SORT_EIGS_BY_ABS_VALUE
        % assert that except for the 3 largest ones, the eigs are sorted by abs value.
        [~,ids] = sort(-abs(evalues(4:end)));
        evalues(4:end) = evalues(3+ids);
        evecs(:,4:end) = evecs(:,3+ids);
    end
end

% Print eigenvalues
log_message( '3NX3N Sync Matrix First %d Eigenvalues:\n%s', n_eigs, num2str(evalues') );
if n_eigs>=4
    log_message('Spectral gap: 3rd_eig/4th_eig = %.2f', evalues(3)/evalues(4));
else
    log_message('Cannot compute spectral gap: Need to use n_eigs>=4');
end
% Cancel symmetrization
% till now we used a symmetrized variant of the weighted Sync matrix,
% thus we didn't get the right eigenvectors. to fix that we just need
% to multiply:
evecs = Dhalf * evecs(:,1:3);

% Normalization
% normalize eigenvectors such that every 3X3 block may be orthogonal
for i = 1:3
    evecs(:,i) = sqrt(N)*evecs(:,i)/norm(evecs(:,i));
end

H = evecs';

% Take Rotations from eigenvectors & test orthogonality
R = zeros(3,3,N);
orthogonality_error = zeros(N,1);
non_orthogonal_counter = 0;
for i = 1:N
    R(:,:,i) = H(:, (3*i-2):(3*i) );
    [U,~,V] = svd(R(:,:,i));
    tmp = U*V.';
    orthogonality_error(i) = norm(R(:,:,i)-tmp,'fro'); % distance from orthogonality
    if (ORTHOGONALITY_THRESHOLD>0) && (orthogonality_error(i) > ORTHOGONALITY_THRESHOLD)
        % count non orthogonal rotation
        non_orthogonal_counter = non_orthogonal_counter + 1;
    end
    % fix rotation to be orthogonal
    R(:,:,i) = tmp;
end
log_message('Non orthogonal rotations (above error threshold of %.2f):\t%d', ORTHOGONALITY_THRESHOLD, non_orthogonal_counter);

end
