%{
INPUT:
    data: data(:,:,i) is ith projection
    c: bandlimit
    R: support radius
    numCompThreads: maximum number of CPU to be used
OUTPUT:
    FBCoeff: FB_coeff(:,i) contains the Fourier-Bessel coefficients for ith image
    B: corresponding Bessel stuff
    n_r: FB parameter

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ FBCoeff , B , n_r ] = getFBCoeff( data , c , R , numCompThreads )

n = size( data , 3 );

n_r = ceil( 4*c*R );

% Fourier-Bessel basis
[ basis , sample_points ] = precomp_fb( n_r , R , c );

% Fourier-Bessel coefficients
data2 = cell( numCompThreads , 1 );
nb = floor( n/numCompThreads );
remain = n - nb*numCompThreads;
for i = 1:remain
    data2{i} = data( : , : , (nb+1)*(i-1)+1:(nb+1)*i );
end
count = (nb+1)*remain;
for i = remain+1:numCompThreads
    data2{i} = data( : , : , count + (i-remain-1)*nb+1: count + (i-remain)*nb );
end
clear data

coeff_pos_k = FBcoeff_nfft( data2 , R , basis , sample_points);

% useful Bessel zeros
load bessel.mat  
B_pos = bessel( bessel(:, 4) <= 2*pi*c*R , : );
clear bessel

% extend to negative k
[FBCoeff,B] = extendFB(coeff_pos_k,B_pos);

end









