% function y=icdft(x)
%
% Aliased inverse discrete Fourier transform (IDFT) of the sequence x.
% The DFT is computed directly using O(n^2) operations.
%
% x    The sequence whose IDFT should be computed. Can be of odd or even length.
%
% Returns the aliased IDFT of the sequence x.
% 
% Yoel Shkolnisky 22/10/01

function y=icdft(x)
m=length(x);
y=zeros(1,m);

for k=lowIdx(m):hiIdx(m)
   acc = 0;
   for j=lowIdx(m):hiIdx(m)
      acc = acc + x(toUnaliasedIdx(j,m))* exp(2*pi*i*j*k/m);
   end;
   y(toUnaliasedIdx(k,m)) = acc/m;
end