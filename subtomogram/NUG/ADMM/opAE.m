%{
INPUT:
    X0,X1: index by (1,1),(2,2),...,(n,n) | (2,1),(3,1),...(n,1) | (3,2),(4,2),...,(n,2) |...
    Xq: quaternion block
    AEq: operator for quaternion equality constraint
    n: number of images
OUTPUT:
    yE0 = AE( X0 )
    yE1 = AE( X1 )
    yEq: quaternion

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ yE0 , yE1 , yEq ] = opAE( X0 , X1 , Xq , AEq , n )

% identity
yE0 = X0(:,1:n);
yE1 = X1(:,1:n);

% quaternion
yEq = AEq*[ Xq ; X0(1,(n+1):end) ; X1(1:4,(n+1):end) ];

end









