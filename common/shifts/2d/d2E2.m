function y=d2E2(deltax,rho,N,idx) 
%
% d2E2(j,deltax,rho,N,idx) compute the Hessian of the function E2. 
%
% Yoel Shkolnisky, January 2014.


ll=fix(N/2);
freqrng=-ll:N-ll-1;
[X,Y]=ndgrid(freqrng,freqrng);

y=zeros(numel(idx),4); % The (symmetric matrix of the) Hessian is 
                       % is stored as a vector containing d2x,
                       % dxdy, dydx, d2y  (in that order).

a=2.*pi./N;
g=exp(1i.*a.*(deltax(1).*X(idx)+deltax(2).*Y(idx)));
y(:,1)=(-(a.*X(idx)).^2).*(g.*conj(E2(deltax,rho,N,idx))+conj(g).*E2(deltax,rho,N,idx))...
    +2*(a.*X(idx)).^2.*ones(size(X(idx)));
y(:,2)=(-(a.^2).*X(idx).*Y(idx)).*(g.*conj(E2(deltax,rho,N,idx))+conj(g).*E2(deltax,rho,N,idx))...
    +2*(a.^2).*X(idx).*Y(idx).*ones(size(X(idx)));
y(:,3)=y(:,2);
y(:,4)=(-(a.*Y(idx)).^2).*(g.*conj(E2(deltax,rho,N,idx))+conj(g).*E2(deltax,rho,N,idx))...
    +2*(a.*Y(idx)).^2.*ones(size(X(idx)));
