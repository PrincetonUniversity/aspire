function [ C ] = clstack2C( clstack,L )
% Store the common line vectors c_ij , i,j=1,...K from the common line
% stack. C(:,i,j)=(c_ij(1), c_ij(2))';
% 
% Lanhui Wang
% Jul 18, 2012
K=size(clstack,1);
% L=max(clstack(:));
fprintf('K=%d, L=%d\n',K,L);

C=zeros(2,K,K);
for k1=1:K;
    k2=(k1+1):K;
    l1 = clstack(k1,k2)-1;
    l2 = clstack(k2,k1)-1;
    l1=l1(:);
    l2=l2(:);
    
    x12 = cos(2*pi*l1/L);
    y12 = sin(2*pi*l1/L);
    x21 = cos(2*pi*l2/L);
    y21 = sin(2*pi*l2/L);
    
    C(1,k1,k2)=x12;
    C(2,k1,k2)=y12;
    C(1,k2,k1)=x21;
    C(2,k2,k1)=y21;
    
end;

end

