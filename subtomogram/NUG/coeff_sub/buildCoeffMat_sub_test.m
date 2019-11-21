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
function coeffMat = buildCoeffMat_sub_test( pf , c, R, bandlmt , t , K_sub, ns, b_sub, n_theta, numCompThreads )

% n = size(pf,3);

% Fourier-Bessel coefficients
[ FBCoeff , B , n_r ] = getFBCoeff( pf , c , R , numCompThreads );
delete(gcp)

% Legendre quadrature
[ r , w ] = lgwt( n_r , 0 , c );

% to projections in Fourier-Bessel space
[ P , ~ ] = FBCoeff_to_image( B , r , c , FBCoeff , bandlmt );

% SO(3) grid...objective function doesn't depend on 'beta'
[ alpha , beta, gamma ] = buildSO3Grid_sub(bandlmt); % switched because we do transform over 'alpha' first

% objective functions on SO(3) grid, checked
%P = reshape(pf,[],n_theta/2*ns,K_sub);
objVal = getObjVal_sub_test( w, P , alpha , beta, gamma , K_sub, ns, b_sub, n_theta, numCompThreads );

% SO(3) Fourier transform
%coeffMat = SO3FT( objVal , bandlmt , t , K_sub );
% build objVal sub
count = 0;
objVal_subI = zeros( 1 , round(K_sub*(K_sub-1)/2) );
objVal_subJ = zeros( 1 , round(K_sub*(K_sub-1)/2) );
for j = 1:(K_sub-1)
for i = (j+1):K_sub
    count = count + 1;
    objVal_subI(count) = i;
    objVal_subJ(count) = j;
end
end

Npair = length(objVal_subI);
objVal = reshape(objVal,2*bandlmt,2*bandlmt,2*bandlmt,Npair);
coeffMat = cell(t+1,1);
for i = 1:t+1
    coeffMat{i,1} = zeros((2*i-1)*length(objVal_subI));
end

beta_weights = zeros(2*bandlmt,1);
for l = 0:(2*bandlmt-1)
    
    val = ( 1 ./ ( 2*(0:(bandlmt-1))+1 ) ) .* sin( pi*(2*l+1) * (2*(0:(bandlmt-1))+1) / (4*bandlmt) );
    val = sum(val);
    val = (1/bandlmt) * sin( pi*(2*l+1) / (4*bandlmt) ) * val;
    
    beta_weights(l+1) = val;
    
end

alpha_2d = pi*(0:2*bandlmt-1)/bandlmt;
[alpha_2d,gamma_2d] = ndgrid(alpha_2d,alpha_2d);
beta_1d = pi*(2*(0:(2*bandlmt-1))+1)/(4*bandlmt);

% fft on alpha and gamma, checked
S2 = zeros(2*t+1,2*t+1,2*bandlmt,Npair);
for k = 0:(2*bandlmt-1)
    for i = 1:Npair   
        S2(:,:,k+1,i) = nufft2d1((2*bandlmt)^2,alpha_2d(:),gamma_2d(:),squeeze(objVal(:,k+1,:,i))...
            ,1,1e-15,2*t+1,2*t+1);  % nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt)
    end
end

% k = 2; i = 3; m = 2; mp = -2;
% s21 = 0;
% for j2 = 1:2*bandlmt
%     fac = exp(1i*mp*(j2-1)*pi/bandlmt);
%     for j1 = 1:2*bandlmt
%         fac1 = exp(1i*m*(j1-1)*pi/bandlmt);
%         s21 = s21+fac1*fac*objVal(j1,k+1,j2,i);
%     end
% end
    
% beta transform and rearrange to output matrix
for l = 0:t
    
    dl = 2*l+1;
    coeffMat_l = zeros(dl*K_sub, dl*K_sub);
    
    % 'beta' transform
    for k = 0:(2*bandlmt-1)
        W_beta = beta_weights(k+1) * (wignerd(l,beta_1d(k+1)))';
        for i = 1:Npair   
            idx_i = (objVal_subI(i)-1)*dl;
            idx_j = (objVal_subJ(i)-1)*dl;
            coeffMat_l(idx_i+1:idx_i+dl,idx_j+1:idx_j+dl) = coeffMat_l(idx_i+1:idx_i+dl,idx_j+1:idx_j+dl)...
                + W_beta.*(S2(t+1-l:t+1+l,t+1-l:t+1+l,k+1,i)');
        
        end
    end 
    
    % ji position is the conjugate transpose
    coeffMat{l+1,1} = coeffMat_l + coeffMat_l';
    
end

% l = 2; m = 1; mp = -1; i = 2;
% f_hat = 0;
% for k = 1:2*bandlmt
%     d = wignerd(l,beta_1d(k));
%     f_hat = f_hat + beta_weights(k)* d(l+1+mp,l+1+m) * S2(t+1+m,t+1+mp,k,i); % why transpose?
% end

end









