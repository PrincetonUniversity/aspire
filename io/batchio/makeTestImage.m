function im=makeTestImage(sz,numkey)
%
% MAKETESTIMAGE     Create test image for testing imagestack
%
% makeTestImage(sz,numkey)
%    Create a image of size sz showing the number in numkey.
%       sz is the size of the image (scalar or a vector for rectangular
%       images). numkey is the number to draw on the image
%       Returns a grayscale image.
%
% Yoel Shkolnisky, May 2016.

if isscalar(sz)
    sz=[sz sz];
end
im=zeros(sz(1),sz(2),3);
text=sprintf('%06d',numkey);
position =  [floor(sz(1)/2) floor(sz(2)/2)];
im = insertText(im,position,text,'AnchorPoint','Center','FontSize',round(min(sz)*0.2));
im=rgb2gray(im);