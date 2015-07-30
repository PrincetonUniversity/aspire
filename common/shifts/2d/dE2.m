function y=dE2(deltax,rho,N,idx) 
%
% dE2(j,deltax,rho,N,idx) compute the gradient of the function E2. 
% The result is an array of size numel(idx) by 2,where the first colums is
% dE2/dx and the second is dE2/dy.
%
% Yoel Shkolnisky, January 2014.

ll=fix(N/2);
freqrng=-ll:N-ll-1;
[X,Y]=ndgrid(freqrng,freqrng);
y=zeros(numel(idx),2);
g=exp(2.*pi.*1i.*(X(idx).*deltax(1)+Y(idx).*deltax(2))./N);
y(:,1)=(g.*conj(E2(deltax,rho,N,idx))-conj(g).*E2(deltax,rho,N,idx)).*2.*pi.*1i.*X(idx)/N;
y(:,2)=(g.*conj(E2(deltax,rho,N,idx))-conj(g).*E2(deltax,rho,N,idx)).*2.*pi.*1i.*Y(idx)/N;