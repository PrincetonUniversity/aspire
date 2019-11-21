%{
INPUT:
    n: number of images
    G0,G1: SO(3) elements corresponding to the inquality operator AI
    C0,C1: coefficient matrix
    yI: dual variable for Fejer kernel
    S0,S1: dual variable for PSD constraint
    Sq: dual variable for quaternion PSD constraint
    AEq: operator for quaternion equality constraint
    AEqAEqtInv: pinv(AE*AE') for quaternion...other are identity
OUTPUT:
    yE0,yE1: dual variable for Xk_ii = identity
    yEq: dual variable for X{1}_ij <-> quaternion

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ yE0 , yE1 , yEq ] = update_yE( n , G0 , G1 , C0 , C1 , yI , S0 , S1 , Sq , AEq , AEqAEqtInv )

[ matI0 , matI1 ] = opAIadj( yI , n , G0 , G1 );

mat0 = C0 - matI0 - S0;
mat1 = C1 - matI1 - S1;

[ yE0 , yE1 , yEq ] = opAE( mat0 , mat1 , -Sq , AEq , n );

yEq = AEqAEqtInv*yEq;

end









