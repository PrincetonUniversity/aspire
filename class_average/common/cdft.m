% function y=cdft(x)
%
% Aliased DFT of the sequence x.
% The DFT is computed directly using O(n^2) operations.
%
% x    The sequence whose DFT should be computed. Can be of odd or even length.
%
% Returns the aliased DFT of the sequence x.
% 
% Yoel Shkolnisky 22/10/01

function y=cdft(x)
m=length(x);
y=zeros(1,m);

for k=lowIdx(m):hiIdx(m)
   acc = 0;
   for j=lowIdx(m):hiIdx(m)
      acc = acc + x(toUnaliasedIdx(j,m))* exp(-2*pi*i*j*k/m);
   end;
   y(toUnaliasedIdx(k,m)) = acc;
end