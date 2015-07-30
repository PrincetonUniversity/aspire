function y=d2E2(deltax,rho,N) 
%
% d2E2(j,deltax,rho,N) compute the Hessian of the function E2. This is
% a two-dimensional version of the function d2E. 
%
% Yoel Shkolnisky, January 2014.

[X,Y]=meshgrid(-N:N,-N:N);
y=zeros(size(X,1),size(X,2),4); % The (symmetric matrix of the) Hessian is 
                                % is stored as a vector containing d2x,
                                % dxdy, dydx, d2y  (in that order).

a=2.*pi./(2*N+1);
g=exp(1i.*a.*(deltax(1).*X+deltax(2).*Y));
y(:,:,1)=(-(a.*X).^2).*(g.*conj(E2(deltax,rho,N))+conj(g).*E2(deltax,rho,N))...
    +2*(a.*X).^2.*ones(size(X));
y(:,:,2)=(-(a.^2).*X.*Y).*(g.*conj(E2(deltax,rho,N))+conj(g).*E2(deltax,rho,N))...
    +2*(a.^2).*X.*Y.*ones(size(X));
y(:,:,3)=y(:,:,2);
y(:,:,4)=(-(a.*Y).^2).*(g.*conj(E2(deltax,rho,N))+conj(g).*E2(deltax,rho,N))...
    +2*(a.*Y).^2.*ones(size(X));
