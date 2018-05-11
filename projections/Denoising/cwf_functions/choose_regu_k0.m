function [regu]=choose_regu_k0(A,w)

% Choose regularization parameter for the least squares system
% w: Weights
% A: cell, CTF operators
% Tejal Bhamre: updated Dec 2015

D=size(A{1});

%disp('here')
for i=1:length(A)
    B{i}=A{i}'*A{i};
    %%E=E+w(i)*A{i}'*C{i}*A{i};
end
[D1,D2]=size(B{1});% D2 dimension of C, D1 dimension of X
AA=zeros(D1^2,D2^2);
for i=1:length(B)
    AA=AA+kron(B{i},B{i})*w(i);  %+kron(conj(A{i}),conj(A{i}))*w(i);
end
ee=eig(AA);
regu=sqrt(min(ee)*max(ee)/2);
