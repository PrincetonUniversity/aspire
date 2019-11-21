function [ mat0_out , mat1_out ] = C_to_CkFull( mat0 , mat1 , d0 , d1 , n , k )

[ Xrow , Xcol ] = Xrow_and_Xcol( n );

mat0_out = zeros( n*d0(k) ); mat0_out = mat0_out(:);
count = 0;
for col = 1:d0(k)
for row = 1:d0(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d0(k),n*d0(k)] , (Xrow-1)*d0(k)+row , (Xcol-1)*d0(k)+col );
    
    mat0_out(idx) = mat0( sum(d0(1:(k-1)).^2) + count , (n+1):end );

end
end
mat0_out = reshape( mat0_out , n*d0(k) , n*d0(k) );
mat0_out = mat0_out + mat0_out';

mat1_out = zeros( n*d1(k) ); mat1_out = mat1_out(:);
count = 0;
for col = 1:d1(k)
for row = 1:d1(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d1(k),n*d1(k)] , (Xrow-1)*d1(k)+row , (Xcol-1)*d1(k)+col );
    
    mat1_out(idx) = mat1( sum(d1(1:(k-1)).^2) + count , (n+1):end );

end
end
mat1_out = reshape( mat1_out , n*d1(k) , n*d1(k) );
mat1_out = mat1_out + mat1_out';

end









