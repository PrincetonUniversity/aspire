%{
INPUT:
    ( X0 , X1 ) = AIadj( yI )
    n: number of signals
    G0,G1: SO(3) elements corresponding to the inquality operator AI
OUTPUT:
    X0,X1: index by (1,1),(2,2),...,(n,n) | (2,1),(3,1),...(n,1) | (3,2),(4,2),...,(n,2) |...

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ X ] = opAIadj_noJ_hetero( yI, yIh , n , G )

% X0 = G0'*yI;
% X0 = [ zeros( size(X0,1) , n ) , X0 ];
% 
% X1 = G1'*yI;
% X1 = [ zeros( size(X1,1) , n ) , X1 ];

X = G'*yI;
X = [zeros(size(X,1),n), X];

X = yIh;

end









