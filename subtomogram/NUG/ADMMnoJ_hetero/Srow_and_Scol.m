%{
INPUT:
    n: number of images
OUTPUT:
    Srow,Scol: useful indices (see <update_S.m>)

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ Srow , Scol ] = Srow_and_Scol( n )

Srow = zeros( round(n*(n-1)/2) , 1 );
Scol = zeros( size(Srow) );

idx1 = 0;
for i = 2:n
    
    idx0 = idx1 + 1;
    idx1 = idx0 + ( n - i );
    
    Srow( idx0:idx1 ) = i:n;
    Scol( idx0:idx1 ) = i-1;
    
end

Srow = [ (1:n)' ; Srow ];
Scol = [ (1:n)' ; Scol ];

end









