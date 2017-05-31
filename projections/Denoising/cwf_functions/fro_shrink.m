function [eigs_shrink]=fro_shrink(eigs, noise_v, gamma)
% Frobenius norm shrinkage to eigenvalues, see Donoho et al. Optimal Shrinkage of Eigenvalues in the Spiked Covariance Model
% http://arxiv.org/abs/1311.0851 Optimal Shrinkage of Eigenvalues in the Spiked Covariance Model
% Tejal Bhamre, Dec 2015
% INPUTS:
%eigs: Eigenvalues
%noise_v: Noise variance
%gamma: p/n
%OUTPUT: Eigenvalues after shrinkage



eigs=diag(eigs);
eigs=eigs/noise_v;
cutoff=(1+sqrt(gamma)).^2;
ll=@(x,gamma) ((x+1-gamma)+sqrt((x+1-gamma).^2-4*x))/2;
cl=@(l,gamma) sqrt((1-gamma*(l-1).^2)./(1+gamma*(l-1))) ;
eta=@(x,cutoff,l,c,gamma) (1 + c.^2.*(l-1)).*(x>cutoff) ;
l=ll(eigs,gamma);
eigs_shrink=eta(eigs,cutoff,l,cl(l,gamma),gamma);
eigs_shrink=noise_v*(max(eigs_shrink-1,0));
end

