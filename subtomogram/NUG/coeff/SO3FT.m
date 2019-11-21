function coeffMat = SO3FT( objVal , bandlmt , t , n )

% 'alpha' and 'gamma' transform
alpha = (pi/bandlmt)*(0:(2*bandlmt-1));
EXP = zeros( 2*t+1 , 2*bandlmt );
for k = (-t):t
    EXP( k+t+1 , : ) = exp( -1i*k*alpha ); % conjugation here
end

SO3FT_objVal = EXP * reshape( objVal , 2*bandlmt , [] );
SO3FT_objVal = reshape( SO3FT_objVal , (2*t+1) , 2*bandlmt , [] );
SO3FT_objVal = permute( SO3FT_objVal , [2,1,3] );
SO3FT_objVal = EXP * reshape( SO3FT_objVal , 2*bandlmt , [] );
SO3FT_objVal = reshape( SO3FT_objVal , 2*t+1 , 2*t+1 , [] );
SO3FT_objVal = 1/((2*bandlmt)^2) * SO3FT_objVal;
% NOTE: no transpose here because integral is over \rho^*...conjugation done earlier

% beta
beta = pi*(2*(0:(2*bandlmt-1))+1)/(4*bandlmt);

beta_weights = zeros(2*bandlmt,1);
for l = 0:(2*bandlmt-1)
    
    val = ( 1 ./ ( 2*(0:(bandlmt-1))+1 ) ) .* sin( pi*(2*l+1) * (2*(0:(bandlmt-1))+1) / (4*bandlmt) );
    val = sum(val);
    val = (1/bandlmt) * sin( pi*(2*l+1) / (4*bandlmt) ) * val;
    
    beta_weights(l+1) = val;
    
end

% useful entries
coeffRow = zeros( round(n*(n-1)/2) , 1 );
coeffCol = zeros( size(coeffRow) );

idx1 = 0;
for i = 2:n
    
    idx0 = idx1 + 1;
    idx1 = idx0 + ( n - i );
    
    coeffRow( idx0:idx1 ) = i:n;
    coeffCol( idx0:idx1 ) = i-1;
    
end

coeffMat = cell(t+1,1);
for k = 0:t
    
    dk = 2*k + 1;
    
    % 'beta' transform
    W_beta = zeros(dk,dk);
    for l = 0:(2*bandlmt-1)
        W_beta = W_beta + beta_weights(l+1) * (wignerd(k,beta(l+1)))';
    end
    
    coeff_k = SO3FT_objVal( ((-k):k)+t+1 , ((-k):k)+t+1 , : ) .* repmat( W_beta , 1 , 1 , round(n*(n-1)/2) );
    coeff_k = reshape( coeff_k , dk^2 , [] );
    
    % re-arrange to matrix format
    coeffMat_k = zeros( n*dk , n*dk );
    coeffMat_k = coeffMat_k(:);
    count = 0;
    for col = 1:dk
    for row = 1:dk
        
        count = count + 1;
        
        idx = sub2ind( [n*dk,n*dk] , (coeffRow-1)*dk+row , (coeffCol-1)*dk+col ); % useful indices
        coeffMat_k(idx) = coeff_k( count , : );
        
    end
    end
    
    coeffMat_k = reshape( coeffMat_k , n*dk , n*dk );
    coeffMat_k = coeffMat_k + coeffMat_k';
    
    coeffMat{k+1} = coeffMat_k;
    
end

end









