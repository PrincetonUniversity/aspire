function y=dE(j,deltax,rho,N) 
%
% dE(j,deltax,rho,N) Derivative of the function E(j,deltax,rho,N).
%
% Yoel Shkolnisky, January 2014

g=exp(2.*pi.*1i.*j.*deltax./(2*N+1));
y=(g.*conj(E(j,deltax,rho,N))-conj(g).*E(j,deltax,rho,N)).*2.*pi.*1i.*j/(2*N+1);
