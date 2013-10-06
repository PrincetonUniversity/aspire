% function y=icfftd(x,d)
%
% 1-D Aliased inverse FFT of the multi-dimensional image x along dimension d.
%
% x   The image whose FFT should be computed. Can be of odd or even length.
% d   The dimension to transform into the frequency domain.
%
% Returns the 1-D aliased inverse FFT of the image x along dimension d.
% 
% Yoel Shkolnisky 13/1/03

function y=icfftd(x,d)
numDims = ndims(x);
s = size(x);
y=zeros(s); % prepare an output array with the same size as the original array.

s(d) = 1; % ignore the size of dimension d

% vecNum is the number of 1-D vectors we need to process. Its value is the 
% product of the elements of s without the dimension d (This is why we set 
% the d'th element to 1).
% For example, if the input array is 4x3x5 and we want to compute inverse fft along the
% third dimension, then we have to compute 4x3 ffts of vectors of length 5.
vecNum = prod(s);

for k=1:vecNum
   % generate the index of the vector to process
   [subscript{1:numDims}] = ind2sub(s,k); 
   
   % set the d'th dimension to ":" so we extract the entire vector along dimension d
   subscript{d} = ':'; 
   
   % compute the Fourier transform of the extracted vector and return it back to its place
   y(subscript{:}) = icfft(x(subscript{:}));
end;