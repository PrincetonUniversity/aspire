function [A_sdp, b_sdp, c_sdp, KK, pars] = SDPLR_prep(K, S)
% prepare sedumi input format
%
% Lanhui Wang, Aug 8, 2013
n = 2*K;
KK.s = n; % semidefinite matrix of size 2K x 2K

c_sdp = reshape(-S,n^2,1); % -S because SDPLR solves a minimization problem

n_eqs = 2*K + K;
%     nnz = n_eqs;
A_sdp=sparse(1:n_eqs,...
    [((1:n)-1)*n+(1:n) ((1:K)-1)*n+K+(1:K)],ones(n_eqs,1),n_eqs,n*n);

b_sdp = [ones(n,1) ; zeros(K,1)];
pars.limit=360000;