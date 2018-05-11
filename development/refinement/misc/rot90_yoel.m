function B=rot90(A)
%
% Rotate the image A by 90 degrees CCW.
%   B=rot90(A);
%
% Yoel Shkolnisky, February 2011.


% In matrix notation
% R90 = [ cosd(90) -sind(90) ;
%         sind(90)  cosd(90) ]
% and then B(p)=A(R90*p) where p=(x,y).'.
% Specifically, B(x,y)=A(-y,x);

B=A.';
B=flipdim(B,1);