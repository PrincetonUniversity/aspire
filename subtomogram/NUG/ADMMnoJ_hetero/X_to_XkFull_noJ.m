function [ mat_out ] = X_to_XkFull_noJ( mat, d , n , k )

[ Xrow , Xcol ] = Xrow_and_Xcol( n );

mat_out = zeros( n*d(k) ); mat_out = mat_out(:);
count = 0;
for col = 1:d(k)
for row = 1:d(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d(k),n*d(k)] , (Xrow-1)*d(k)+row , (Xcol-1)*d(k)+col );
    
    mat_out(idx) = mat( sum(d(1:(k-1)).^2) + count , (n+1):end );

end
end
mat_out = reshape( mat_out , n*d(k) , n*d(k) );
mat_out = mat_out + mat_out' + eye(size(mat_out));

% mat1_out = zeros( n*d1(k) ); mat1_out = mat1_out(:);
% count = 0;
% for col = 1:d1(k)
% for row = 1:d1(k)
% 
%     count = count + 1;
% 
%     % useful indices
%     idx = sub2ind( [n*d1(k),n*d1(k)] , (Xrow-1)*d1(k)+row , (Xcol-1)*d1(k)+col );
%     
%     mat1_out(idx) = mat1( sum(d1(1:(k-1)).^2) + count , (n+1):end );
% 
% end
% end
% mat1_out = reshape( mat1_out , n*d1(k) , n*d1(k) );
% mat1_out = mat1_out + mat1_out' + eye(size(mat1_out));

end









