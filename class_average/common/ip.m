% function c=ip(a,b);
%
% Inner product of 2 multi-dimensional arrays.
% For two matrices the inner product is given by
%
%         N1   N2
%    c = sum (sum a(k,l)*conj(b(k,l))
%        k=1  l=1 
%
% Yoel Shkolnisky 9/2/02

function c=ip(a,b);
	c = sum(conj(a(:)).*b(:));