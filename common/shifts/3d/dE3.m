function y=dE3(deltax,rho,N,idx) 
%
% dE3(j,deltax,rho,N,idx) compute the gradient of the function E3.
% The result is an array of size numel(idx) by 3,where the first colums is
% dE3/dx, the second is dE3/dy, and the third is dE3/dz.
%
% Yoel Shkolnisky, January 2014.

ll=fix(N/2);
freqrng=-ll:N-ll-1;
[X,Y,Z]=ndgrid(freqrng,freqrng,freqrng);
y=zeros(numel(idx),3);
g=exp(2.*pi.*1i.*(X(idx).*deltax(1)+Y(idx).*deltax(2)+Z(idx).*deltax(3))./N);
y(:,1)=(g.*conj(E3(deltax,rho,N,idx))-conj(g).*E3(deltax,rho,N,idx)).*2.*pi.*1i.*X(idx)/N;
y(:,2)=(g.*conj(E3(deltax,rho,N,idx))-conj(g).*E3(deltax,rho,N,idx)).*2.*pi.*1i.*Y(idx)/N;
y(:,3)=(g.*conj(E3(deltax,rho,N,idx))-conj(g).*E3(deltax,rho,N,idx)).*2.*pi.*1i.*Z(idx)/N;