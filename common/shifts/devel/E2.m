function y=E2(deltax,rho,N) 
%
% E2(j,deltax,rho,N) two-dimensional version of the function E.
%
% Yoel Shkolnisky, January 2014.

[X,Y]=meshgrid(-N:N,-N:N);
y=(exp(2*pi*1i*(X.*deltax(1)+Y.*deltax(2))/(2*N+1))-rho);