function y=dE2(deltax,rho,N) 
%
% dE2(j,deltax,rho,N) compute the gradient of the function E2. This is
% a two-dimensional version of the function dE. 
%
% Yoel Shkolnisky, January 2014.

[X,Y]=meshgrid(-N:N,-N:N);
y=zeros(size(X,1),size(X,2),2);
g=exp(2.*pi.*1i.*(X.*deltax(1)+Y.*deltax(2))./(2*N+1));
y(:,:,1)=(g.*conj(E2(deltax,rho,N))-conj(g).*E2(deltax,rho,N)).*2.*pi.*1i.*X/(2*N+1);
y(:,:,2)=(g.*conj(E2(deltax,rho,N))-conj(g).*E2(deltax,rho,N)).*2.*pi.*1i.*Y/(2*N+1);