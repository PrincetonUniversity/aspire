% function y=cfftn(x)
%
% Aliased n-dimensional FFT of the image x.
% The FFT is computed using O((n^d)logn) operations, where d is the dimension of the image.
%
% x   The image whose FFT should be computed. Can be of odd or even length in each dimension.
%
% Returns the aliased n-D FFT of the image x.
% 
% Yoel Shkolnisky 11/1/03


function y=cfftn(x)
    y = fftshift(fftn(ifftshift(x)));
