%{
INPUT:
    G0,G1: SO(3) elements corresponding to the inquality operator AI
OUTPUT:
    lambda: largest eigenvalue of AI*AI'

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function lambda = findLambda( G0 , G1 )

lambda = eigs(G0*G0'+G1*G1',1);

end









