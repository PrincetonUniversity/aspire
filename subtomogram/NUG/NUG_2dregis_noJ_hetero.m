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
%function [ X_out , runtime ] = NUG_2dregis_noJ_hetero( numIter , coeff )

if length(coeff) < 2
    error('NUG_cryoEM: length(coeff) < 2')
end

n = size(coeff{1},1);
t = length(coeff) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fejer kernel

tS2 = 20; % spherical t-design (may be extremal) for Hopf fibration
[ d, G, bI ] = buildG_Fejer_noJ( t , tS2 );
d = [1;d];
bI = repmat( bI , size(G,1) , round(n*(n-1)/2) );

if exist(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'_noJ.mat'],'file') == 2
    load(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'_noJ.mat'])
else
    lambda = real(findLambda_noJ(G));
    save(['NUG/lambda/lambda_t',num2str(t),'_Fejer_Hopf',num2str(tS2),'_noJ.mat'],'lambda')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficient matrix (without 0th representation)

C = cell(t+1,1);

Cnorm = 0; Xnorm = 0;

for k = 0:t
    
    dk = 2*k+1;
    
    [ T , Tinv ] = realY_to_complexY(k);
    
    C{k+1} = real( kron(eye(n),Tinv) * coeff{k+1} * kron(eye(n),T) );
    
    Cnorm = Cnorm + norm(C{k+1},'fro')^2;
    Xnorm = Xnorm + dk*n^2;
    
end
%clear('coeff')

Cnorm = sqrt( Cnorm );
Xnorm = sqrt( Xnorm );

% normalize
for k = 0:t
    C{k+1} = ( Xnorm / Cnorm ) * C{k+1};
end

% stack C{k}
C_tmp = zeros( sum(d.^2), round(n*(n+1)/2));
for k = 1:t+1
    
    [ C_k ] = CkFull_to_Ck_noJ( C{k} , d , n , k ); % C_k's columns are C_k(i,j)
                                                % vectorized    
    C_tmp( (sum(d(1:(k-1)).^2)+1):sum(d(1:k).^2) , : ) = C_k;
    
end
%clear('C')

% split into C0 and C1
% [ d0 , d1 , G0 , G1 , bI ] = buildG_Fejer( t , tS2 );
% C0 = zeros( sum(d0.^2) , round(n*(n+1)/2) );
% C1 = zeros( sum(d1.^2) , round(n*(n+1)/2) );
% for k = 1:t
%     
%     [ C0_k , C1_k ] = splitBlocks( C{k} , k , n );
%     
%     [ C0_k , C1_k ] = CkFull_to_Ck( C0_k , C1_k , d0 , d1 , n , k );
%     
%     C0( (sum(d0(1:(k-1)).^2)+1):sum(d0(1:k).^2) , : ) = C0_k;
%     C1( (sum(d1(1:(k-1)).^2)+1):sum(d1(1:k).^2) , : ) = C1_k;
%     
% end
%clear('C')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADMM

fprintf('starting ADMM...\n')

tic;

load('NUG/ADMMnoJ/AE_noJ.mat')
AEq = AE_noJ;
AEqAEqtInv = AEAEtInv_noJ ;

% initialize change X0 and X1
[ X , Xq ] = initializeX_noJ( n , d );
[ S, Sq ] = initializeS_noJ( n , d );
yI = zeros( size(G,1) , round(n*(n-1)/2) );
[ yE, yEq ] = update_yE_noJ( n , G, C_tmp, yI , S, Sq , AEq , AEqAEqtInv );

rho = 1;
count_rhoRepeat = 0;

for iter = 1:numIter
    
    fprintf(['iteration #',num2str(iter),'\n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADMM blocks
    
    [ S , Sq ] = update_S_noJ( n , d , rho , G , C_tmp , yE , yEq , yI , X , Xq , AEq );
    
    [ yE , yEq ] = update_yE_noJ( n , G , C_tmp , yI , S , Sq , AEq , AEqAEqtInv );
    
    yI = update_yI_withSemiProximal_noJ( n , rho , lambda , bI , G , C_tmp , yE , yEq , yI , S , X , AEq );
    
    [ yE , yEq ] = update_yE_noJ( n , G , C_tmp , yI , S , Sq , AEq , AEqAEqtInv );
    
    [ X , Xq ] = update_X_noJ( n , rho , G , C_tmp , yE , yEq , yI , S , Sq , X , Xq , AEq );
    
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
% X = cell(t,1);
% for k = 1:t
%     [ X0k , X1k ] = X_to_XkFull( X0 , X1 , d0 , d1 , n , k );
%     X{k} = combineBlocks( X0k , X1k , k , n );
% end

% change X into cell structure
X_out = cell(t,1);
for k = 1:t
    [X_out{k}] = X_to_XkFull_noJ(X, d, n, k);
end

runtime = toc;

fprintf('DONE!\n')

%end









