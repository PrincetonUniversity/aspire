function [R, eigenvalues, orthogonality_error] = cryo_sync3n_S_to_rot (S, W, Dhalf, n_eigs, ORTHOGONALITY_THRESHOLD)
% CRYO_SYNC3N_S_TO_ROT  Extract rotations from synchronization matrix.
%
% cryo_sync3n_S_to_rot (S, n_eigs, ORTHOGONALITY_THRESHOLD)
%   Given a 3Nx3N synchronization matrix S, whose (i,j)'th block of size
%   3x3 is an estimate Ru'*Rj, for estimate the rortaions Ri.
%   The function also compute the top n_eigs eigenvalues of S, which can be
%   used to inspect the spectral gap.
%
% cryo_sync3n_S_to_rot (S, W, Dhalf, N, n_eigs, ORTHOGONALITY_THRESHOLD)
%   Weight each block in S according to the weights in W.
%
% Input:
%   S   3NX3N Synchronization matrix whose blocks are the relative
%       rotations, as returned from cryo_sync3n_syncmatrix.
%   W   3NX3N weights matrix (without rows normalization).
%   Dhalf   Normalizing matrix of W
%   n_eigs  Number of eigenvalues to compute
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
%          ||Ri-orthogonalize(Ri)||_F.
%
%XXX Wny W is not NxN by rather 3Nx3N?
% Why do we need both W and Dhalf? Why not computing D from W?
%
% This function was revised from
%   rotations (S, W, Dhalf, N, n_eigs, ORTHOGONALITY_THRESHOLD)
% Parameter N was removed and is calculated from the dimensions of S.
%
% Yoel Shkolnisky, March 2016.
% Based on code by Ido Greenberg, 2015

if nargin<=3
    if isscalar(W)
        % The parameters are (S, n_eigs, ORTHOGONALITY_THRESHOLD)
        if nargin==3
            ORTHOGONALITY_THRESHOLD=Dhalf;
        end
        n_eigs=W;
        W=1;
        Dhalf=1;
    end
elseif nargin<5
    ORTHOGONALITY_THRESHOLD=-1; % Ignore
end

% Input validation
if n_eigs < 3
    error('n_eigs must be >= 3 integer! otherwise H cannot be deduced from Sync Matrix. Currently n_eigs = %d', n_eigs);
end

% Compute eigenvectors
[evecs, evalues] = eigs( Dhalf * (W.*S) * Dhalf , n_eigs);

% eigs() returns eigenvalues sorted by absolute value. we wish to sort
% the first 3 by real value, and the rest by absolute value.
% GGG actually it's not trivial that we prefer abs upon real value
evalues = diag(evalues);
[~,ids] = sort(-evalues); % (assert 3 largest eigenvalues are first)
evalues = evalues(ids);
evecs = evecs(:,ids);
if n_eigs > 3
    [~,ids] = sort(-abs(evalues(4:end))); % (assert other eigs are sorted by abs)
    evalues(4:end) = evalues(3+ids);
    evecs(:,4:end) = evecs(:,3+ids);
end
eigenvalues = evalues;

% Print eigenvalues
log_message( '3NX3N Sync Matrix First %d Eigenvalues:\n%s\n', n_eigs, num2str(evalues') );
if n_eigs>=4
    log_message('Spectral gap: 3rd_eig/4th_eig = %.2f\n', evalues(3)/evalues(4));
else
    log_message('Cannot compute spectral gap. Used n_eigs>=4');
end
% Cancel symmetrization
% till now we used a symmetrized variant of the weighted Sync matrix,
% thus we didn't get the right eigenvectors. to fix that we just need
% to multiply:
evecs = Dhalf * evecs(:,1:3);

% Normalization
% normalize eigenvectors such that every 3X3 block may be orthogonal
N=size(S,1)/3;
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
    if (orthogonality_error(i) > ORTHOGONALITY_THRESHOLD) && (ORTHOGONALITY_THRESHOLD>0)
        % count non orthogonal rotation
        non_orthogonal_counter = non_orthogonal_counter + 1;
    end
    % fix rotation to be orthogonal
    R(:,:,i) = tmp;
end
log_message('Non orthogonal rotations (above error threshold of %.2f):\t%d\n', ORTHOGONALITY_THRESHOLD, non_orthogonal_counter);

end
