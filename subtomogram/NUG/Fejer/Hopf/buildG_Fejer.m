function [ d0 , d1 , G0 , G1 , bI ] = buildG_Fejer( t , t_S2 )

SO3 = discretizeSO3(t_S2);
nG = length(SO3);

d0 = zeros(t,1);
d1 = zeros(t,1);
for k = 1:t
    d0(k) = k;
    d1(k) = k+1;
end
G0 = zeros( nG , sum(d0.^2) );
G1 = zeros( nG , sum(d1.^2) );

for m = 1:nG
    
    alpha = SO3(m,1); beta = SO3(m,2); gamma = SO3(m,3);
    
    row0_m = zeros(1,sum(d0.^2));
    row1_m = zeros(1,sum(d1.^2));
    
    for k = 1:t
        
        tmp = (1/2)*(t-k+2)*(t-k+1)*(k+1/2)*realWignerD(k,alpha,beta,gamma);
        
        tmp0 = reshape( tmp(2:2:(2*k),2:2:(2*k)) , (d0(k))^2 , 1 );
        row0_m( (sum(d0(1:(k-1)).^2)+1):(sum(d0(1:k).^2)) ) = tmp0';
        
        tmp1 = reshape( tmp(1:2:(2*k+1),1:2:(2*k+1)) , (d1(k))^2 , 1 );
        row1_m( (sum(d1(1:(k-1)).^2)+1):(sum(d1(1:k).^2)) ) = tmp1';
        
    end
    
    G0(m,:) = row0_m;
    G1(m,:) = row1_m;
    
end

bI = ( -(1/2)*(t+2)*(t+1)*(1/2) );

end









