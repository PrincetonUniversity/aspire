function B=rot180(A)
%
% Rotate the image A by 180 degrees CCW.
%   B=rot180(A);
%
% Yoel Shkolnisky, February 2011.

 
% In matrix notation
% R180 = [ cosd(180) -sind(180) ;
%          sind(180)  cosd(180) ]
% and then B(p)=A(R180*p) where p=(x,y).'.
% Specifically, B(x,y)=A(-x,-y);

B=flipdim(A,1);
B=flipdim(B,2);