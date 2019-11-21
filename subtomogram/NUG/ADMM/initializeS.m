%{
INPUT:
    n: number of images
    d0,d1: vector containing dimensions of the representations
OUTPUT:
    S0: k-by-k block
    S1: (k+1)-by-(k+1) block
    Sq: quaternion block

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ S0 , S1 , Sq ] = initializeS( n , d0 , d1 )

S0 = zeros(sum(d0.^2),round(n*(n+1)/2));
S1 = zeros(sum(d1.^2),round(n*(n+1)/2));
Sq = zeros(16,round(n*(n-1)/2));

end









