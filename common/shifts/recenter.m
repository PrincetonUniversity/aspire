function [cX,cm]=recenter(X)
% RECENTER  Recenter an image or volume according to its center of mass.
%
% [cX,cm]=recenter(X)
%   Recenter the image/volume X according to its center of mass. X must be
%   a square (image) or a cube (volume). Returns the centered image/volume
%   cX as well as the center of mass of the input X.
%
% Yoel Shkolnisky, November 2014.

% Check input
dims=ndims(X);
if dims~=2 && dims~=3
    error('X must be an square image or a cube volume');
end
if size(X,1)~=size(X,2)
    error('All dimensions of X must be the same');
end
if (dims==3) && size(X,1)~=size(X,3)
   error('All dimensions of X must be the same');
end 

cm=CenterOfMass(X);
if dims==2
    cX=reshift_image(X,cm);
else
    cX=reshift_vol(X,cm);
end

% Compute center of mass again to make sure the object is centered.
% If it is not, it is an indication for a problem.
cm2=CenterOfMass(cX);
if norm(cm2)>0.1
    warning('Could not center object');
end
