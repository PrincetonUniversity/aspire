function [est_inv_rots, Gram] = est_orientations_LUD(common_lines_matrix, L, pars)
% Determine orientations of rotations using LUD and semidefinite relaxation
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
%   pars: options to use different solvers and the spectral norm constraint.
%       The spectral norm constraint
%       ||G||_2 <= alpha*K is used if pars.alpha is in the range
%       [2/3, 1). Default: pars.alpha = 0, no spectral norm constraint is
%       used. pars.solver = 'ADMM', the ADMM solver is used.
%
% Output:
%   est_inv_rots:  a 3 x 3 x K matrix, the ith estimated inverse rotation
%       matrix is est_inv_rots(:,:,i).
%
%   Gram: The 2K x 2K Gram matrix before the rounding procedure.
%
% Lanhui Wang, Aug 8, 2013

if (nargin < 3) || (~isfield(pars, 'alpha'));    pars.alpha = 0;   end
if (nargin < 3) || (~isfield(pars, 'solver'));    pars.solver = 'ADMM'; end
% Compute the 2K x 2K commonline matrix S
K = size(common_lines_matrix, 1);

if max(strcmp(pars.solver, {'IRLS', 'ADMM'})) < 1
    error('%s is not valid! Please choose to use ADMM or IRLS',pars.solver);

elseif strcmp(pars.solver, 'ADMM')
    fprintf('Starting ADMM\n');
    C = clstack2C( common_lines_matrix,L );
    opts.tol = 1e-3;
    if (pars.alpha < 2/3) || (pars.alpha >= 1)
        fprintf('No spectral norm constraint is used.\n');
        Gram = cryoEMSDPL12N_vsimple(K, eye(2*K), C,opts);
    else
        fprintf('Spectral norm constraint is used.\n');
        alpha = pars.alpha;
        Gram = cryoEMSDPL12N(K, eye(2*K), alpha*K, C,opts);
    end
    Gram=[Gram(1:2:2*K,1:2:2*K) Gram(1:2:2*K,2:2:2*K);...
    Gram(2:2:2*K,1:2:2*K) Gram(2:2:2*K,2:2:2*K)];

elseif strcmp(pars.solver, 'IRLS')
    fprintf('Starting IRLS\n');
    if ~isfield(pars,'NumOfIters'); pars.NumOfIters = 10; end % default number
    % of iterations is 10.
    N = pars.NumOfIters;
    
    if ~isfield(pars,'epsilon'); pars.epsilon = 1e-3; end % default epsilon
    epsilon = pars.epsilon;
    
    S = construct_S(common_lines_matrix, L);
    W = ones(2*K);
    Gram = eye(2*K);
    if (pars.alpha < 2/3) || (pars.alpha >= 1)
        fprintf('No spectral norm constraint is used.\n');
        [A, b, c, KK, pars] = SDPLR_prep(K, S);
        x_sdp = [];
        y = [];
        info = [];
        r = [];
        n = 2*K;
        for k = 1 : N
            [x_sdp,~,~,~] = sdplr(A,b,W(:).*c,KK,pars,[],x_sdp,y,info,r);
            Gram = reshape(x_sdp(1:n^2),n,n);
            [W, res] = update_weights(S, Gram, K, epsilon);
        end
    else
        fprintf('Spectral norm constraint is used.\n');
        opt.ctyep=0;
        alpha = pars.alpha;
        for k = 1 : N
            [Gram, ~, ~, ~, ~] = cryoEMADM(K, Gram, alpha*K, W.*S,opt);
            [W, res] = update_weights(S, Gram, K, epsilon);
        end
    end
end
% Using the rounding procedure to recover the rotations from the
% Gram matrix.
fprintf('Startig the rounding procedure\n');
est_inv_rots = deterministic_rounding(Gram);
end

function [W, res] = update_weights(S, Gram, K, epsilon)

W=S.*Gram;
weights=W(1:K,1:K)+W(1:K,K+1:2*K)+W(K+1:2*K,1:K)+W(K+1:2*K,K+1:2*K);
weights=sqrt(abs(2-2*weights));
res=sum(weights(:));
weights=1./sqrt(weights.^2+epsilon^2);
W=[weights weights; weights weights];

end

