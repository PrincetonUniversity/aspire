function [ P , theta ] = FBCoeff_to_image( B , r , c , FBCoeff , bandlmt )

n = size( FBCoeff , 2 );

% Bessel basis
Bessel_basis0 = zeros( length(B(:,3)) , length(r) );
for m = 1:length(r)
    Bessel_basis0(:,m) = besselj( B(:,1) , B(:,3)*r(m)/c ) ./ ( c * sqrt(pi) * abs( besselj( B(:,1)+1 , B(:,3) ) ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_FB = size(Bessel_basis0,1);
n_r = size(Bessel_basis0,2);

theta = (pi/bandlmt)*(0:(2*bandlmt-1)); % row vector...same as 'alpha' and 'gamma'
n_theta = length(theta);

% Fourier
k = reshape( B(:,1) , 1 , 1 , n_FB );
k = repmat( k , n_r , n_theta , 1 );

Fourier_basis = repmat( theta , n_r , 1 , n_FB );
Fourier_basis = k.*Fourier_basis;
Fourier_basis = exp( 1i*Fourier_basis );

% Bessel
Bessel_basis = reshape( Bessel_basis0' , n_r , 1 , n_FB );
Bessel_basis = repmat( Bessel_basis , 1 , n_theta , 1 );

% Fourier-Bessel
FB_basis = Fourier_basis.*Bessel_basis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate at ( r , theta )

P = zeros( n_r , n_theta , n );

for i = 1:n
    
    P_i = reshape( FBCoeff(:,i) , 1 , 1 , length(FBCoeff(:,i)) );
    P_i = repmat( P_i , n_r , n_theta , 1 );
    
    P(:,:,i) = sum( P_i.*FB_basis , 3 );
    
end

end









