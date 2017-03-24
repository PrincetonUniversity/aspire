function [est_inv_rots, Gram] = est_orientations_LS(common_lines_matrix, L, pars)
% Determine orientations of rotations using LS and semidefinite relaxation
%
% Input:
%   common_lines_matrix: a K x K matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2.
%
%   L:  the number of radial lines within the Fourier transform of a
%       projection.
%
%   pars: option to use the spectral norm constraint. The spectral norm
%       constraint ||G||_2 <= alpha*K is used if pars.alpha is in the range
%       [2/3, 1). Default: pars.alpha = 0, no spectral norm constraint is
%       used.
%
% Output:
%   est_inv_rots:  a 3 x 3 x K matrix, the ith estimated inverse rotation
%       matrix is est_inv_rots(:,:,i).
%
%   Gram: The 2K x 2K Gram matrix before the rounding procedure.
%
% Lanhui Wang, Aug 8, 2013

if (nargin < 3) || (~isfield(pars, 'alpha'));    pars.alpha = 0;   end
% Compute the 2K x 2K commonline matrix S
K = size(common_lines_matrix, 1);
S = construct_S(common_lines_matrix, L);
if (pars.alpha < 2/3) || (pars.alpha >= 1)
    fprintf('No spectral norm constraint is used.\n Startig SDPLR');
    [A, b, c, KK, pars] = SDPLR_prep(K, S);
    [x_sdp,~,~,~] = sdplr(A,b,c,KK,pars);
    n = 2*K;
    Gram = reshape(x_sdp(1:n^2),n,n);
else
    fprintf('Spectral norm constraint is used.\n Startig ADMM');
    opt.ctyep=0;
    alpha = pars.alpha;
    [Gram, ~, ~, ~, ~] = cryoEMADM(K, eye(2*K), alpha*K, S,opt);
end

% Using the rounding procedure to recover the rotations from the
% Gram matrix.
fprintf('Startig the rounding procedure');
%est_inv_rots = randomized_rounding(Gram);
est_inv_rots = deterministic_rounding(Gram);