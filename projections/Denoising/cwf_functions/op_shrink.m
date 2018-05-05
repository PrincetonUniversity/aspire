% OP_SHRINK Apply operator norm shrinkage
%
% Usage
%    eigs_shrink = op_shrink(eigs, noise_v, gamma);
%
% Input
%    eigs: The eigenvalues in diagonal matrix form.
%    noise_v: The noise variance.
%    gamma: The aspect ratio (p/n) of the data matrix.
%
% Output
%    eigs_shrink: The eigenvalues after shrinkage.
%
% Description
%    For details, seee Donoho et al., "Optimal Shrinkage of Eigenvalues in the
%    Spiked Covariance Model".

% Written by Tejal Bhamre - Dec 2015
% Reformatted and documented by Joakim Anden - 2018-Apr-19

function eigs_shrink = op_shrink(eigs, noise_v, gamma)
    eigs=diag(eigs);
    eigs=eigs/noise_v;

    cutoff=(1+sqrt(gamma)).^2;

    ll=@(x,gamma)(((x+1-gamma)+sqrt((x+1-gamma).^2-4*x))/2);

    l=ll(eigs,gamma);

    eta=@(x,cutoff,l)(l.*(x>cutoff));
    eigs_shrink=eta(eigs,cutoff,l);
    eigs_shrink=noise_v*(max(eigs_shrink-1,0));
end
