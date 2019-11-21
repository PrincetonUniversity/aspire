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
function coeffMat = buildCoeffMat( projs , c , R , bandlmt , t , numCompThreads )

n = size(projs,3);

% Fourier-Bessel coefficients
[ FBCoeff , B , n_r ] = getFBCoeff( projs , c , R , numCompThreads );
delete(gcp)

% Legendre quadrature
[ r , w ] = lgwt( n_r , 0 , c );

% to projections in Fourier-Bessel space
[ P , theta ] = FBCoeff_to_image( B , r , c , FBCoeff , bandlmt );

% SO(3) grid...objective function doesn't depend on 'beta'
[ gamma , alpha ] = buildSO3Grid(bandlmt); % switched because we do transform over 'alpha' first

% objective functions on SO(3) grid
objVal = getObjVal( w , P , theta , alpha , gamma , numCompThreads );

% SO(3) Fourier transform
coeffMat = SO3FT( objVal , bandlmt , t , n );

end









