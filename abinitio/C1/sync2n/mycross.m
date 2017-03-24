function c=mycross(a,b)
%
% a and b must be two column vectors.
%
% Yoel Shkolnisky, August 2010.

c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
 
