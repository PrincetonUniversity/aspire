function [ mat_out ] = S_to_SkFull_noJ( mat , d , n , k )

% mat0(:,1:n) = 0.5*mat0(:,1:n);
% mat1(:,1:n) = 0.5*mat1(:,1:n);
% 
% [ Srow , Scol ] = Srow_and_Scol( n );
% 
% %
% mat0_out = zeros( n*d0(k) ); mat0_out = mat0_out(:);
% count = 0;
% for col = 1:d0(k)
% for row = 1:d0(k)
% 
%     count = count + 1;
% 
%     % useful indices
%     idx = sub2ind( [n*d0(k),n*d0(k)] , (Srow-1)*d0(k)+row , (Scol-1)*d0(k)+col );
%     
%     mat0_out(idx) = mat0( sum(d0(1:(k-1)).^2) + count , : );
% 
% end
% end
% mat0_out = reshape( mat0_out , n*d0(k) , n*d0(k) );
% mat0_out = mat0_out + mat0_out';
% 
% %
% mat1_out = zeros( n*d1(k) ); mat1_out = mat1_out(:);
% count = 0;
% for col = 1:d1(k)
% for row = 1:d1(k)
% 
%     count = count + 1;
% 
%     % useful indices
%     idx = sub2ind( [n*d1(k),n*d1(k)] , (Srow-1)*d1(k)+row , (Scol-1)*d1(k)+col );
%     
%     mat1_out(idx) = mat1( sum(d1(1:(k-1)).^2) + count , : );
% 
% end
% end
% mat1_out = reshape( mat1_out , n*d1(k) , n*d1(k) );
% mat1_out = mat1_out + mat1_out';

% add for no J ambiguity
mat(:,1:n) = 0.5*mat(:,1:n);

[ Srow , Scol ] = Srow_and_Scol( n );

%
mat_out = zeros( n*d(k) ); mat_out = mat_out(:);
count = 0;
for col = 1:d(k)
for row = 1:d(k)

    count = count + 1;

    % useful indices
    idx = sub2ind( [n*d(k),n*d(k)] , (Srow-1)*d(k)+row , (Scol-1)*d(k)+col );
    
    mat_out(idx) = mat( sum(d(1:(k-1)).^2) + count , : );

end
end
mat_out = reshape( mat_out , n*d(k) , n*d(k) );
mat_out = mat_out + mat_out';

end









