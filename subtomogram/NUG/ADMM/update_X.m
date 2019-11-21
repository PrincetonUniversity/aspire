%{
INPUT:
    n: number of images
    rho: ADMM parameter
    G0,G1: SO(3) elements corresponding to the inquality operator AI
    C0,C1: coefficient matrix
    yE0,yE1: dual variable for Xk_ii = identity
    yEq: dual variable for X{1}_ij <-> quaternion
    yI: dual variable for Fejer kernel
    S0,S1: dual variable for PSD constraint
    Sq: dual variable for quaternion PSD constraint
    X0_in,X1_in: primal variable
    Xq_in: primal variable for quaternion
    AEq: operator for quaternion equality constraint
OUTPUT:
    X0_out,X1_out: primal variable
    Xq_out: primal variable for quaternion

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ X0_out , X1_out , Xq_out ] = update_X( n , rho , G0 , G1 , C0 , C1 , yE0 , yE1 , yEq , yI , S0 , S1 , Sq , X0_in , X1_in , Xq_in , AEq )

[ matE0 , matE1 , matEq ] = opAEadj( yE0 , yE1 , yEq , AEq , n );

[ matI0 , matI1 ] = opAIadj( yI , n , G0 , G1 );

X0_out = X0_in + rho*( matE0 + matI0 + S0 - C0 );
X1_out = X1_in + rho*( matE1 + matI1 + S1 - C1 );
Xq_out = Xq_in + rho*( matEq + Sq );

end









