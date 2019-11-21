%{
INPUT:
    numIter: number of ADMM iterations
    coeff: Fourier coefficients of the objective function
OUTPUT:
    X: SDP solution
    runtime: runtime

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ X , runtime ] = NUG_cryoEM( numIter , coeff )

if length(coeff) < 2
    error('NUG_cryoEM: length(coeff) < 2')
end

n = size(coeff{1},1);
t = length(coeff) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fejer kernel

tS2 = 20; % spherical t-design (may be extremal) for Hopf fibration
[ d0 , d1 , G0 , G1 , bI ] = buildG_Fejer( t , tS2 );
bI = repmat( bI , size(G0,1) , round(n*(n-1)/2) );

if exist(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'.mat'],'file') == 2
    load(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'.mat'])
else
    lambda = real(findLambda(G0,G1));
    save(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'.mat'],'lambda')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficient matrix (without 0th representation)

C = cell(t,1);

Cnorm = 0; Xnorm = 0;

for k = 1:t
    
    dk = 2*k+1;
    
    [ T , Tinv ] = realY_to_complexY(k);
    
    C{k} = real( kron(eye(n),Tinv) * coeff{k+1} * kron(eye(n),T) );
    
    Cnorm = Cnorm + norm(C{k},'fro')^2;
    Xnorm = Xnorm + dk*n^2;
    
end
clear('coeff')

Cnorm = sqrt( Cnorm );
Xnorm = sqrt( Xnorm );

% normalize
for k = 1:t
    C{k} = ( Xnorm / Cnorm ) * C{k};
end

% split into C0 and C1
C0 = zeros( sum(d0.^2) , round(n*(n+1)/2) );
C1 = zeros( sum(d1.^2) , round(n*(n+1)/2) );
for k = 1:t
    
    [ C0_k , C1_k ] = splitBlocks( C{k} , k , n );
    
    [ C0_k , C1_k ] = CkFull_to_Ck( C0_k , C1_k , d0 , d1 , n , k );
    
    C0( (sum(d0(1:(k-1)).^2)+1):sum(d0(1:k).^2) , : ) = C0_k;
    C1( (sum(d1(1:(k-1)).^2)+1):sum(d1(1:k).^2) , : ) = C1_k;
    
end
clear('C')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADMM

fprintf('starting ADMM...\n')

tic;

load('NUG/ADMM/AE.mat')
AEq = AE;
AEqAEqtInv = AEAEtInv ;

% initialize
[ X0 , X1 , Xq ] = initializeX( n , d0 , d1 );
[ S0 , S1 , Sq ] = initializeS( n , d0 , d1 );
yI = zeros( size(G0,1) , round(n*(n-1)/2) );
[ yE0 , yE1 , yEq ] = update_yE( n , G0 , G1 , C0 , C1 , yI , S0 , S1 , Sq , AEq , AEqAEqtInv );

rho = 1;
count_rhoRepeat = 0;

for iter = 1:numIter
    
    fprintf(['iteration #',num2str(iter),'\n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADMM blocks
    
    [ S0 , S1 , Sq ] = update_S( n , d0 , d1 , rho , G0 , G1 , C0 , C1 , yE0 , yE1 , yEq , yI , X0 , X1 , Xq , AEq );
    
    [ yE0 , yE1 , yEq ] = update_yE( n , G0 , G1 , C0 , C1 , yI , S0 , S1 , Sq , AEq , AEqAEqtInv );
    
    yI = update_yI_withSemiProximal( n , rho , lambda , bI , G0 , G1 , C0 , C1 , yE0 , yE1 , yEq , yI , S0 , S1 , X0 , X1 , AEq );
    
    [ yE0 , yE1 , yEq ] = update_yE( n , G0 , G1 , C0 , C1 , yI , S0 , S1 , Sq , AEq , AEqAEqtInv );
    
    [ X0 , X1 , Xq ] = update_X( n , rho , G0 , G1 , C0 , C1 , yE0 , yE1 , yEq , yI , S0 , S1 , Sq , X0 , X1 , Xq , AEq );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 'rho' update
    
    if iter < 50
        rho = rho/1.01;
    else
        count_rhoRepeat = count_rhoRepeat + 1;
        if count_rhoRepeat > 5
            count_rhoRepeat = 0;
            rho = rho/1.01;
        end
    end
    
end

% combine X0 and X1
X = cell(t,1);
for k = 1:t
    [ X0k , X1k ] = X_to_XkFull( X0 , X1 , d0 , d1 , n , k );
    X{k} = combineBlocks( X0k , X1k , k , n );
end

runtime = toc;

fprintf('DONE!\n')

end









