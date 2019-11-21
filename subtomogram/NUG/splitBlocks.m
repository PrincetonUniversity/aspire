%{
INPUT:
    mat: (2k+1)-by-(2k+1) block
    k: degree of representation
    n: number of images
OUTPUT:
    mat0: k-by-k block
    mat1: (k+1)-by-(k+1) block

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ mat0 , mat1 ] = splitBlocks( mat , k , n )

dk = 2*k+1;

mat0 = zeros(k*n,k*n);
mat1 = zeros((k+1)*n,(k+1)*n);

for i = 1:n
for j = 1:n
    
    mat0(((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)) = mat(((i-1)*dk+2):2:(i*dk-1),((j-1)*dk+2):2:(j*dk-1));
    mat1(((i-1)*(k+1)+1):(i*(k+1)),((j-1)*(k+1)+1):(j*(k+1))) = mat(((i-1)*dk+1):2:(i*dk),((j-1)*dk+1):2:(j*dk));
    
end
end

end









