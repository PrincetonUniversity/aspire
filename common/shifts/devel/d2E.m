function y=d2E(j,deltax,rho,N) 
% d2E(j,deltax,rho,N) second derivative of E(j,deltax,rho,N).
%
% Yoel Shkolnisky, January 2014.

a=2.*pi.*j/(2*N+1);
g=exp(1i.*a.*deltax);
y=(-a.^2).*(g.*conj(E(j,deltax,rho,N))+conj(g).*E(j,deltax,rho,N))...
    +2*a.^2.*ones(size(j));
