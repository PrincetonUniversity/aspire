function y=E(j,deltax,rho,N) 
%
% E(j,deltax,rho,N) error in phases induced by the estimated shifts deltax
% Take as an input the estimated shift deltax and the measured phase
% difference between the signals, and check how well the phases induced by
% deltax agree with the measured phases rho. That is, return the difference
% between the induced phases and the measured ones.
%
% Yoel Shkolnisky, January 2014.

y=(exp(2*pi*1i*j*deltax/(2*N+1))-rho);