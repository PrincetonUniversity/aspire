%{
INPUT:
    w: weights for numerical integration
    P: projections in Fourier space
    theta: angle discretization
    alpha: 1st Euler angle
    gamma: 3rd Euler angle
OUTPUT:
    objVal: value of the objective function at (alpha,gamma)X(i,j)

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function objVal = getObjVal( w , P , theta , alpha , gamma , numCompThreads )

n = size(P,3);

% theta_i = gamma - pi/2 and theta_j = -alpha - pi/2......THIS HAS BEEN VERIFIED
[ ~ , idx_I ] = min( transpose( abs( repmat(theta,length(gamma),1) - mod( repmat(gamma,1,length(theta)) - pi/2 , 2*pi ) ) ) );
[ ~ , idx_J ] = min( transpose( abs( repmat(theta,length(alpha),1) - mod( -repmat(alpha,1,length(theta)) - pi/2 , 2*pi ) ) ) );

% build objVal sub
count = 0;
objVal_subI = zeros( 1 , round(n*(n-1)/2) );
objVal_subJ = zeros( 1 , round(n*(n-1)/2) );
for j = 1:(n-1)
for i = (j+1):n
    count = count + 1;
    objVal_subI(count) = i;
    objVal_subJ(count) = j;
end
end

% get objVal entries corresponding to sub
objVal = zeros( length(idx_I) , length(objVal_subI) );

parpool('local',numCompThreads);
parfor m = 1:length(idx_I)
    
    % L1
    objVal_ij = abs( P( : , idx_I(m) , objVal_subI ) - P( : , idx_J(m) , objVal_subJ ) );
    objVal_ij = sum( objVal_ij .* repmat( w , 1 , 1 , size(objVal_ij,3) ) , 1 );
    
    objVal( m , : ) = squeeze(objVal_ij);
    
end
delete(gcp);

end









