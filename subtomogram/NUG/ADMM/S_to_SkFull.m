function [ mat0_out , mat1_out ] = S_to_SkFull( mat0 , mat1 , d0 , d1 , n , k )

mat0(:,1:n) = 0.5*mat0(:,1:n);
mat1(:,1:n) = 0.5*mat1(:,1:n);

[ Srow , Scol ] = Srow_and_Scol( n );

%
mat0_out = zeros( n*d0(k) ); mat0_out = mat0_out(:);
count = 0;
for col = 1:d0(k)
for row = 1:d0(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d0(k),n*d0(k)] , (Srow-1)*d0(k)+row , (Scol-1)*d0(k)+col );
    
    mat0_out(idx) = mat0( sum(d0(1:(k-1)).^2) + count , : );

end
end
mat0_out = reshape( mat0_out , n*d0(k) , n*d0(k) );
mat0_out = mat0_out + mat0_out';

%
mat1_out = zeros( n*d1(k) ); mat1_out = mat1_out(:);
count = 0;
for col = 1:d1(k)
for row = 1:d1(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d1(k),n*d1(k)] , (Srow-1)*d1(k)+row , (Scol-1)*d1(k)+col );
    
    mat1_out(idx) = mat1( sum(d1(1:(k-1)).^2) + count , : );

end
end
mat1_out = reshape( mat1_out , n*d1(k) , n*d1(k) );
mat1_out = mat1_out + mat1_out';

end









