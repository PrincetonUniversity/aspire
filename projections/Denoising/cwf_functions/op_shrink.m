
function [eigs_shrink]=op_shrink(eigs, noise_v, gamma)
% Operator norm shrinkage to eigenvalues, see Donoho et al. Optimal Shrinkage of Eigenvalues in the Spiked Covariance Model
% Tejal Dec 2015
% INPUTS:
%eigs: Eigenvalues
%noise_v: Noise variance
%gamma: p/n
%OUTPUT: Eigenvalues after shrinkage

eigs=diag(eigs);
eigs=eigs/noise_v;
cutoff=(1+sqrt(gamma)).^2;
ll=@(x,gamma) ((x+1-gamma)+sqrt((x+1-gamma).^2-4*x))/2;
l=ll(eigs,gamma);
eta=@(x,cutoff,l) l.*(x>cutoff) ;
eigs_shrink=eta(eigs,cutoff,l);
eigs_shrink=noise_v*(max(eigs_shrink-1,0));
end

