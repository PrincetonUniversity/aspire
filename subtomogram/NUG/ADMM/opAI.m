%{
INPUT:
    X0,X1: index by (1,1),(2,2),...,(n,n) | (2,1),(3,1),...(n,1) | (3,2),(4,2),...,(n,2) |...
    n: number of signals
    G0,G1: SO(3) elements corresponding to the inquality operator AI
OUTPUT:
    yI = AI( X0 , X1 )

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function yI = opAI( X0 , X1 , n , G0 , G1 )

yI = G0*X0(:,(n+1):end) + G1*X1(:,(n+1):end);

end









