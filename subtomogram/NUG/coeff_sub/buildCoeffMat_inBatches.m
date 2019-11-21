%{
INPUT:
    data: data(:,:,i) is ith projection
    varNoise: variance of the noise
    bandlmt: bandlimit of the objective function
    t: truncate at (t+1)th degree representation
    numCompThreads: maximum number of CPU to be used
OUTPUT:
    coeffMat: Fourier coefficients of the objective function

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function coeffMat = buildCoeffMat_inBatches( projs , c , R , bandlmt , t , numCompThreads , numBatch )

n = size(projs,3);

% Fourier-Bessel coefficients
[ FBCoeff , B , n_r ] = getFBCoeff( projs , c , R , numCompThreads );
delete(gcp)

% Legendre quadrature
[ r , w ] = lgwt( n_r , 0 , c );

% to projections in Fourier-Bessel space
[ P , theta ] = FBCoeff_to_image( B , r , c , FBCoeff , bandlmt );

% SO(3) grid...objective function doesn't depend on 'beta'
alpha = (pi/bandlmt)*(0:(2*bandlmt-1));
alpha = alpha';
alpha = repmat( alpha , 2*bandlmt , 1 );

gamma = (pi/bandlmt)*(0:(2*bandlmt-1));
gamma = gamma';
gamma = kron( gamma , ones(2*bandlmt,1) );

% theta_i = gamma - pi/2 and theta_j = -alpha - pi/2...THIS HAS BEEN VERIFIED
[ ~ , idx_I ] = min(transpose(abs( repmat(theta,length(gamma),1) - mod( repmat(gamma,1,length(theta)) - pi/2 , 2*pi ))));
[ ~ , idx_J ] = min(transpose(abs( repmat(theta,length(alpha),1) - mod( -repmat(alpha,1,length(theta)) - pi/2 , 2*pi ))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize coefficient matrix
coeffMat = cell(t+1,1);
for k = 0:t
    dk = 2*k + 1;
    coeffMat{k+1} = zeros(n*dk,n*dk);
end

% objVal sub
count = 0;
objVal_subI = zeros( round(n*(n-1)/2) , 1 );
objVal_subJ = zeros( round(n*(n-1)/2) , 1 );
for j = 1:(n-1)
for i = (j+1):n
    count = count + 1;
    objVal_subI(count) = i;
    objVal_subJ(count) = j;
end
end

sizeBatch = floor(length(objVal_subI)/numBatch);
if sizeBatch == 0
    error('buildCoeffMat_inBatches: sizeBatch is 0...not enough images')
end

idx_start = sizeBatch*((1:(numBatch-1))-1) + 1;
idx_end = sizeBatch*(1:(numBatch-1));

% iterate over columns of coefficient matrix
fprintf('building SO(3) coefficients in batches...\n')
for b = 1:(length(idx_start)+1)
    
    fprintf(['batch #',num2str(b),'\n'])
    
    if b == (length(idx_start)+1) % last batch
        
        % objective functions on SO(3) grid
        objVal = getObjVal_inBatches( w , P , idx_I , idx_J , objVal_subI((idx_end(b-1)+1):end) , objVal_subJ((idx_end(b-1)+1):end) );

        % SO(3) Fourier transform
        coeffMat_b = SO3FT_inBatches( objVal , bandlmt , t , n , objVal_subI((idx_end(b-1)+1):end) , objVal_subJ((idx_end(b-1)+1):end) );
        
        clear objVal
        
    else
        
        % objective functions on SO(3) grid
        objVal = getObjVal_inBatches( w , P , idx_I , idx_J , objVal_subI(idx_start(b):idx_end(b)) , objVal_subJ(idx_start(b):idx_end(b)) );

        % SO(3) Fourier transform
        coeffMat_b = SO3FT_inBatches( objVal , bandlmt , t , n , objVal_subI(idx_start(b):idx_end(b)) , objVal_subJ(idx_start(b):idx_end(b)) );
        
        clear objVal
        
    end
    
    % update coeffMat
    for k = 0:t
        coeffMat{k+1} = coeffMat{k+1} + coeffMat_b{k+1};
    end
    
    clear coeffMat_b

end
fprintf('DONE!\n')

end









