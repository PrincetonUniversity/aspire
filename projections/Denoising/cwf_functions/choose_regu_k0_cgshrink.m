function [regu]=choose_regu_k0_cgshrink(A,inv_sqrt_Epop,w)
% Choosing the regularization parameter for the k=0 block while solving CG

D=size(A{1});

%disp('here')
for i=1:length(A)
B1{i}=inv_sqrt_Epop*A{i}'*A{i};
B2{i}=A{i}'*A{i}*inv_sqrt_Epop;
B{i}=A{i}'*A{i};

end

[D1,D2]=size(B1{1});% D2 dimension of C, D1 dimension of X
AA=zeros(D1^2,D2^2);
AA1=zeros(D1^2,D2^2);

for i=1:length(B1)
AA=AA+kron(B2{i}',B1{i})*w(i);  %+kron(conj(A{i}),conj(A{i}))*w(i);
AA1=AA1+kron(B{i}',B{i})*w(i);  %+kron(conj(A{i}),conj(A{i}))*w(i);
end
sprintf('Condition number is %f', cond(AA))
ee=eig(AA);
regu=sqrt(min(ee)*max(ee)/2);
