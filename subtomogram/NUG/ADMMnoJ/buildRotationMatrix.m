function rotationMatrix = buildRotationMatrix( X_out )

dk = 3;

Xk = X_out{1};

n = round( size(Xk,1) / dk );

rotationMatrix = cell(n,n);

[ T , Tinv ] = realY_to_complexY( 1 );

for i = 2:n
for j = 1:(i-1)
    
    Rij = Xk( ((i-1)*dk+1):(i*dk) , ((j-1)*dk+1):(j*dk) );
    
    Wij = T * Rij * Tinv;
    
    % get 'alpha' and 'gamma' up to PI
    gamma_ij = ( angle( Wij(1,dk) ) + angle( Wij(dk,dk) ) ) / 2;
    alpha_ij = ( angle( Wij(dk,1) ) + angle( Wij(dk,dk) ) ) / 2;
    beta_ij = acos(Wij((dk+1)/2,(dk+1)/2));
    
%     theta_i = -alpha_ij - pi/2; theta_j = gamma_ij - pi/2; % this is the transpose...possibly the convention in synchronization?
%     
%     idx_i = n_theta * theta_i / (2*pi);
%     idx_j = n_theta * theta_j / (2*pi);
    
    rotationMatrix{i,j} = EulerZYZ_to_rot(alpha_ij,beta_ij,gamma_ij);
    rotationMatrix{j,i} = transpose(rotationMatrix{i,j});
    
end
end

end









