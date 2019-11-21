%{
INPUT:
    mat0: k-by-k block
    mat1: (k+1)-by-(k+1) block
    k: degree of representation
    n: number of images
OUTPUT:
    mat: (2k+1)-by-(2k+1) block

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function mat = combineBlocks( mat0 , mat1 , k , n )

dk = 2*k+1;

mat = zeros(dk*n,dk*n);

for i = 1:n
for j = 1:n
    
    mat(((i-1)*dk+2):2:(i*dk-1),((j-1)*dk+2):2:(j*dk-1)) = mat0(((i-1)*k+1):(i*k),((j-1)*k+1):(j*k));
    mat(((i-1)*dk+1):2:(i*dk),((j-1)*dk+1):2:(j*dk)) = mat1(((i-1)*(k+1)+1):(i*(k+1)),((j-1)*(k+1)+1):(j*(k+1)));
    
end
end

end









