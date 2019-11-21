%{
INPUT:
    yE0 = AE( X0 )
    yE1 = AE( X1 )
    yEq: quaternion
    AEq: adjoint of AEQ (see <opAE.m>)
    n: number of images
OUTPUT:
    X0: X0{k} = X0k
    X1: X1{k} = X1k
    Xq: quaternion block

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ X0 , X1 , Xq ] = opAEadj( yE0 , yE1 , yEq , AEq , n )

X0 = zeros( size(yE0,1) , round(n*(n+1)/2) );
X1 = zeros( size(yE1,1) , round(n*(n+1)/2) );

X0(:,1:n) = yE0;
X1(:,1:n) = yE1;

tmp = AEq'*yEq;
Xq = tmp( 1:16 , : );
X0(1,(n+1):end) = tmp( 17 , : );
X1(1:4,(n+1):end) = tmp( 18:21 , : );

end









