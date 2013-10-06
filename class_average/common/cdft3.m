% function y=cdft3(x)
%
% Aliased 3D DFT of the image x.
% The DFT is computed directly using O(n^6) operations.
% The function serves as a reference function to test the correctness of the 3D
% Fourier routines.
%
% x   The sequence whose DFT should be computed. Can be of odd or even length.
%
% Returns the aliased 3D DFT of the image x.
% 
% Yoel Shkolnisky 11/1/03

function y=cdft3(x)
m=size(x);
y=zeros(m);

for xi1=lowIdx(m(1)):hiIdx(m(1))       %omegaX direction
   for xi2=lowIdx(m(2)):hiIdx(m(2))    %omegaY direction
      for xi3=lowIdx(m(3)):hiIdx(m(3)) %omegaZ direction
         acc = 0;
         for u=lowIdx(m(1)):hiIdx(m(1))       % x direction
            for v=lowIdx(m(2)):hiIdx(m(2))    % y direction
               for w=lowIdx(m(3)):hiIdx(m(3)) % z direction
                  im_coord = toUnaliasedCoord([u,v,w],m);               
                  acc = acc + x(im_coord{:})* exp(-2*pi*i*(xi1*u/m(1)+xi2*v/m(2)+xi3*w/m(3)));
               end           
             end      
          end
          freq_coord = toUnaliasedCoord([xi1,xi2,xi3],m);
          y(freq_coord{:}) = acc;
       end 
   end
end
