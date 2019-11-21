%{
INPUT:
    n: number of images
    d0,d1: vector containing dimensions of the representations
    rho: ADMM parameter
    G0,G1: SO(3) elements corresponding to the inquality operator AI
    C0,C1: coefficient matrix
    yE0,yE1: dual variable for Xk_ii = identity
    yEq: dual variable for X{1}_ij <-> quaternion
    yI: dual variable for Fejer kernel
    X0,X1: primal variable
    Xq: primal variable for quaternion
    AEq: operator for quaternion equality constraint
OUTPUT:
    S0,S1: dual variable for PSD constraint
    Sq: dual variable for quaternion PSD constraint

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ S0 , S1 , Sq ] = update_S_parallelEig( n , d0 , d1 , rho , G0 , G1 , C0 , C1 , yE0 , yE1 , yEq , yI , X0 , X1 , Xq , AEq , eig_numThreads )

[ matE0 , matE1 , matEq ] = opAEadj( yE0 , yE1 , yEq , AEq , n );

[ matI0 , matI1 ] = opAIadj( yI , n , G0 , G1 );

%
S0 = C0 - matE0 - matI0 - (1/rho)*X0;
S1 = C1 - matE1 - matI1 - (1/rho)*X1;

clear('matE0','matE1','matI0','matI1')

% expand to full matrix form
S0Full = cell(length(d0),1);
S1Full = cell(length(d0),1);
for k = 1:length(d0)
    [ S0_k , S1_k ] = S_to_SkFull( S0 , S1 , d0 , d1 , n , k );
    S0Full{k} = S0_k;
    S1Full{k} = S1_k;
end

% eigen-decomposition
if length(eig_numThreads) > length(d0)
    eig_numThreads = eig_numThreads((end-length(d0)+1):end);
elseif length(eig_numThreads) < length(d0)
    eig_numThreads = [ ones(length(d0)-length(eig_numThreads),1) ; eig_numThreads ];
end

parfor k = 1:length(d0)
    
    warning('off','MATLAB:maxNumCompThreads:Deprecated');
    maxNumCompThreads( eig_numThreads(k) );
    
    [U,D] = eig(S0Full{k});
    S0Full{k} = real( U*subplus(real(D))*U' );
    
    [U,D] = eig(S1Full{k});
    S1Full{k} = real( U*subplus(real(D))*U' );
    
end

maxNumCompThreads(72);

% compress from full matrix form
for k = 1:length(d0)
    [ S0_k , S1_k ] = SkFull_to_Sk( S0Full{k} , S1Full{k} , d0 , d1 , n , k );
    S0( (sum(d0(1:(k-1)).^2)+1):sum(d0(1:k).^2) , : ) = S0_k;
    S1( (sum(d1(1:(k-1)).^2)+1):sum(d1(1:k).^2) , : ) = S1_k;
end

%
Sq = zeros(size(matEq));

for i = 1:size(Xq,2)
    
    [U,D] = eig( reshape( -matEq(:,i) - (1/rho)*Xq(:,i) , 4 , 4 ) );
    
    Sq(:,i) = reshape( real( U*subplus(real(D))*U' ) , 16 , 1 );
    
end

end









