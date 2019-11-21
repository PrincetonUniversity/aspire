%{
INPUT:
    bandlmt: bandlimit
OUTPUT:
    alpha,gamma: grid for SO(3) (objective funcation is constant on 'beta')

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ alpha , gamma ] = buildSO3Grid( bandlmt )

alpha = (pi/bandlmt)*(0:(2*bandlmt-1)); alpha = alpha';
gamma = (pi/bandlmt)*(0:(2*bandlmt-1)); gamma = gamma';

alpha = kron(alpha,ones(2*bandlmt,1));
gamma = repmat(gamma,2*bandlmt,1);

end









