function [ mat_out ] = CkFull_to_Ck_noJ( mat , d , n , k )

[ Xrow , Xcol ] = Xrow_and_Xcol( n );

mat = mat(:);

mat_out = zeros( d(k)^2 , round(n*(n-1)/2) );
count = 0;
for col = 1:d(k)
for row = 1:d(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d(k),n*d(k)] , (Xrow-1)*d(k)+row , (Xcol-1)*d(k)+col );
    
    mat_out(count,:) = mat(idx)';

end
end
mat_out = [ zeros( d(k)^2 , n ) , mat_out ];

% mat0 = mat0(:);
% mat1 = mat1(:);
% 
% mat0_out = zeros( d0(k)^2 , round(n*(n-1)/2) );
% count = 0;
% for col = 1:d0(k)
% for row = 1:d0(k)
% 
%     count = count + 1;
% 
%     % useful indices
%     idx = sub2ind( [n*d0(k),n*d0(k)] , (Xrow-1)*d0(k)+row , (Xcol-1)*d0(k)+col );
%     
%     mat0_out(count,:) = mat0(idx)';
% 
% end
% end
% mat0_out = [ zeros( d0(k)^2 , n ) , mat0_out ];
% 
% mat1_out = zeros( d1(k)^2 , round(n*(n-1)/2) );
% count = 0;
% for col = 1:d1(k)
% for row = 1:d1(k)
% 
%     count = count + 1;
% 
%     % useful indices
%     idx = sub2ind( [n*d1(k),n*d1(k)] , (Xrow-1)*d1(k)+row , (Xcol-1)*d1(k)+col );
%     
%     mat1_out(count,:) = mat1(idx)';
% 
% end
% end
% mat1_out = [ zeros( d1(k)^2 , n ) , mat1_out ];

end









