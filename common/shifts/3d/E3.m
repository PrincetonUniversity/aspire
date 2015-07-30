function y=E3(deltax,rho,N,idx) 
%
% E2(j,deltax,rho,N,idx) difference between the phases induced by the estimated
% shifts deltax and the measured phases rho.
%
% Take as an input the estimated shift deltax and the measured phase
% difference between the signals, and check how well the phases induced by
% deltax agree with the measured phases rho. That is, return the difference
% between the induced phases and the measured ones.
%
% Only the phases whose index is given by idx are considered.
% The sum of the squared absolute value of E3, that is,
%   sum(abs(E3(x,rhat(idx),(n-1)./2,idx)).^2)
% is the L2 error in the phases. This is the expression we minimize (over
% deltax) to find the relative shift between the images.
%
% Yoel Shkolnisky, January 2014.

ll=fix(N/2);
freqrng=-ll:N-ll-1;
[X,Y,Z]=ndgrid(freqrng,freqrng,freqrng);

y=(exp(2*pi*1i*(X(idx).*deltax(1)+Y(idx).*deltax(2)+Z(idx).*deltax(3))/N)-rho);