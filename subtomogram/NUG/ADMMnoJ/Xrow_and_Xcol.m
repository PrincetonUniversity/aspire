%{
INPUT:
    n: number of images
OUTPUT:
    Xrow,Xcol: useful indices

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ Xrow , Xcol ] = Xrow_and_Xcol( n )

Xrow = zeros( round(n*(n-1)/2) , 1 );
Xcol = zeros( size(Xrow) );

idx1 = 0;
for i = 2:n
    
    idx0 = idx1 + 1;
    idx1 = idx0 + ( n - i );
    
    Xrow( idx0:idx1 ) = i:n;
    Xcol( idx0:idx1 ) = i-1;
    
end

end









