function m=RotMatrix2(theta)
% function m=RotMatrix2(theta)
% Compute m such that m*v1 takes the vector v1 into a vector rotated ccw 
% by theta.

c=cos(theta);
s=sin(theta);
m=[c -s; s c];
