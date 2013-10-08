% function y=cdft2(x)
%
% Aliased 2D DFT of the image x.
% The DFT is computed directly using O(n^4) operations.
%
% x   The sequence whose DFT should be computed. Can be of odd or even length.
%
% Returns the aliased 2D DFT of the image x.
% 
% Yoel Shkolnisky 22/10/01

function y=cdft2(x)
m=size(x);
y=zeros(size(x));

for xi1=lowIdx(m(2)):hiIdx(m(2))    %omegaX direction
   for xi2=lowIdx(m(1)):hiIdx(m(1)) %omegaY direction
	   acc = 0;
      for u=lowIdx(m(2)):hiIdx(m(2))    % x direction
         for v=lowIdx(m(1)):hiIdx(m(1)) % y direction
            acc = acc + x(toUnaliasedIdx(v,m(1)),toUnaliasedIdx(u,m(2)))* exp(-2*pi*i*(u*xi1/m(2)+v*xi2/m(1)));
         end      
   	end
      y(toUnaliasedIdx(xi2,m(1)),toUnaliasedIdx(xi1,m(2))) = acc;
   end
end
