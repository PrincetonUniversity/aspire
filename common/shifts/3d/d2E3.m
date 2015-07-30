function y=d2E3(deltax,rho,N,idx) 
%
% d2E3(j,deltax,rho,N,idx) compute the Hessian of the function E3. 
%
% Yoel Shkolnisky, January 2014.

ll=fix(N/2);
freqrng=-ll:N-ll-1;
[X,Y,Z]=ndgrid(freqrng,freqrng,freqrng);

y=zeros(numel(idx),9); %  The (symmetric matrix of the) Hessian is 
                       % is stored as a vector containing 
                       % d2x,  dxdy, dxdz 
                       % dydx, d2y,  dydz
                       % dzdx, dzdy, d2z

a=2.*pi./N;
g=exp(1i.*a.*(deltax(1).*X(idx)+deltax(2).*Y(idx)+deltax(3).*Z(idx)));
y(:,1)=(-(a.*X(idx)).^2).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.*X(idx)).^2.*ones(size(X(idx)));
y(:,2)=(-(a.^2).*X(idx).*Y(idx)).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.^2).*X(idx).*Y(idx).*ones(size(X(idx)));
y(:,3)=(-(a.^2).*X(idx).*Z(idx)).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.^2).*X(idx).*Z(idx).*ones(size(X(idx)));
y(:,4)=y(:,2);
y(:,5)=(-(a.*Y(idx)).^2).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.*Y(idx)).^2.*ones(size(X(idx)));
y(:,6)=(-(a.^2).*Y(idx).*Z(idx)).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.^2).*Y(idx).*Z(idx).*ones(size(X(idx)));
y(:,7)=y(:,3);
y(:,8)=y(:,6);
y(:,9)=(-(a.*Z(idx)).^2).*(g.*conj(E3(deltax,rho,N,idx))+conj(g).*E3(deltax,rho,N,idx))...
    +2*(a.*Z(idx)).^2.*ones(size(X(idx)));
