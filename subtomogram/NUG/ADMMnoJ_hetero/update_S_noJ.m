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
function [ S , Sq ] = update_S_noJ( n , d , rho , G , C , yE , yEq , yI , X , Xq , AEq )

% [ matE0 , matE1 , matEq ] = opAEadj( yE0 , yE1 , yEq , AEq , n );
% 
% [ matI0 , matI1 ] = opAIadj( yI , n , G0 , G1 );
% 
% %
% S0 = C0 - matE0 - matI0 - (1/rho)*X0;
% S1 = C1 - matE1 - matI1 - (1/rho)*X1;
% 
% clear('matE0','matE1','matI0','matI1')
% 
% for k = length(d0):(-1):1
%     
%     [ S0_k , S1_k ] = S_to_SkFull( S0 , S1 , d0 , d1 , n , k );
%     
%     [U,D] = eig(S0_k);
%     S0_k = real( U*subplus(real(D))*U' );
%     
%     [U,D] = eig(S1_k);
%     S1_k = real( U*subplus(real(D))*U' );
%     
%     [ S0_k , S1_k ] = SkFull_to_Sk( S0_k , S1_k , d0 , d1 , n , k );
%     
%     S0( (sum(d0(1:(k-1)).^2)+1):sum(d0(1:k).^2) , : ) = S0_k;
%     S1( (sum(d1(1:(k-1)).^2)+1):sum(d1(1:k).^2) , : ) = S1_k;
%     
%     clear('U','D','S0_k','S1_k')
%     
% end
% 
% %
% Sq = zeros(size(matEq));
% 
% for i = 1:size(Xq,2)
%     
%     [U,D] = eig( reshape( -matEq(:,i) - (1/rho)*Xq(:,i) , 4 , 4 ) );
%     
%     Sq(:,i) = reshape( real( U*subplus(real(D))*U' ) , 16 , 1 );
%     
% end

% add for no J ambiguity
[ matE , matEq ] = opAEadj_noJ( yE , yEq , AEq , n );

[ matI ] = opAIadj_noJ( yI , n , G );

%
S = C - matE - matI - (1/rho)*X;

clear('matE','matI')

for k = length(d):(-1):1
    
    [ S_k ] = S_to_SkFull_noJ( S , d , n , k );
    
    [U,D] = eig(S_k);
    S_k = real( U*subplus(real(D))*U' ); % project onto PSD
    
    [ S_k ] = SkFull_to_Sk_noJ( S_k , d , n , k );
    
    S( (sum(d(1:(k-1)).^2)+1):sum(d(1:k).^2) , : ) = S_k;
    
    clear('U','D','S_k')
    
end

%
Sq = zeros(size(matEq));

for i = 1:size(Xq,2)
    
    [U,D] = eig( reshape( -matEq(:,i) - (1/rho)*Xq(:,i) , 4 , 4 ) );
    
    Sq(:,i) = reshape( real( U*subplus(real(D))*U' ) , 16 , 1 );
    
end

end









