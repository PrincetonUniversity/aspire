function [ mat0_out , mat1_out ] = SkFull_to_Sk( mat0 , mat1 , d0 , d1 , n , k )

mat0 = mat0(:);
mat1 = mat1(:);

[ Srow , Scol ] = Srow_and_Scol( n );

%
mat0_out = zeros( d0(k)^2 , round(n*(n+1)/2) );
count = 0;
for col = 1:d0(k)
for row = 1:d0(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d0(k),n*d0(k)] , (Srow-1)*d0(k)+row , (Scol-1)*d0(k)+col );
    
    mat0_out(count,:) = mat0(idx)';

end
end

%
mat1_out = zeros( d1(k)^2 , round(n*(n+1)/2) );
count = 0;
for col = 1:d1(k)
for row = 1:d1(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d1(k),n*d1(k)] , (Srow-1)*d1(k)+row , (Scol-1)*d1(k)+col );
    
    mat1_out(count,:) = mat1(idx)';

end
end

end









