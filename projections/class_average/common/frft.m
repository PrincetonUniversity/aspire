% function y=frft(x,alpha)
%
% Aliased fractional Fourier transform (FRFT) of the sequence x.
% The FRFT is computed directly using O(n^2) operations.
%
% x       The sequence whose FRFT should be computed. Can be of odd or even length.
% alpha	  The alpha parameter of the fractional Fourier transform
%
% Returns the aliased FRFT with parameter alpha of the sequence x.
% The result is always a row vector.
%
% The fractional Fourier transform y of the sequence x (with parameter alpha) is defined by
%                   n/2-1
%       y(k) =       sum  x(u)*exp(-2*pi*i*k*u*alpha/N),  -n/2 <= k <= n/2-1, N=length(x).
%                   u=-n/2
% The value of the fractional Fourier transform (y) for index k (-n/2 <= k <= n/2-1) is stored in 
% the array y in index toUnalisedIdx(k,N) (which is between 1 and N).
% 
% Yoel Shkolnisky 22/10/01

function y=frft(x,alpha)
m=length(x);
y=zeros(1,m);

for k=lowIdx(m):hiIdx(m)
   acc = 0;
   for j=lowIdx(m):hiIdx(m)
      acc = acc + x(toUnaliasedIdx(j,m))* exp(-2*pi*i*j*k*alpha/m);
   end;
   y(toUnaliasedIdx(k,m)) = acc;
end