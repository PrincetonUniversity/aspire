function B=rot270(A)
%
% Rotate the image A by 270 degrees CCW.
%   B=rot270(A);
%
% Yoel Shkolnisky, February 2011.

 
% In matrix notation
% R270 = [ cosd(270) -sind(270) ;
%          sind(270)  cosd(270) ]
% and then B(p)=A(270*p) where p=(x,y).'.
% Specifically, B(x,y)=A(y,-x);

B=A.';
B=flipdim(B,2);