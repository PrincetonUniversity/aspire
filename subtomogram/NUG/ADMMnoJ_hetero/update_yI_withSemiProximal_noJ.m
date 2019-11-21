%{
INPUT:
    n: number of images
    rho: ADMM parameter
    lambda: eigenvalue for semi-proximal norm
    bI: inequality constants
    G0,G1: SO(3) elements corresponding to the inquality operator AI
    C0,C1: coefficient matrix
    yE0,yE1: dual variable for Xk_ii = identity
    yEq: dual variable for X{1}_ij <-> quaternion
    yI_in: dual variable for Fejer kernel
    S0,S1: dual variable for PSD constraint
    X0,X1: primal variable
    AEq: operator for quaternion equality constraint
OUTPUT:
    yI_out: dual variable for Fejer kernel

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function yI_out = update_yI_withSemiProximal_noJ( n , rho , lambda , bI , G , C , yE , yEq , yI_in , S , X , AEq )

% [ matE0 , matE1 , ~ ] = opAEadj( yE0 , yE1 , yEq , AEq , n );
% 
% [ matI0 , matI1 ] = opAIadj( yI_in , n , G0 , G1 );
% 
% mat0 = C0 - matE0 - matI0 - S0 - (1/rho)*X0;
% mat1 = C1 - matE1 - matI1 - S1 - (1/rho)*X1;
% 
% yI = yI_in + (1/(rho*lambda))*bI + (1/lambda)*opAI( mat0 , mat1 , n , G0 , G1 );
% 
% % nonnegative projection
% yI_out = subplus( yI );

% add for no J ambiguity
[ matE , ~ ] = opAEadj_noJ( yE , yEq , AEq , n );

[ matI ] = opAIadj_noJ( yI_in , n , G );

mat = C - matE - matI - S - (1/rho)*X;

yI = yI_in + (1/(rho*lambda))*bI + (1/lambda)*opAI_noJ( mat , n , G );

% nonnegative projection
yI_out = subplus( yI );

end









