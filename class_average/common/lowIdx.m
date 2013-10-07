% function idx=lowIdx(n)
%
% Return the minimal index for an aliased sequence of length n.
%
% 	n	 the length of the indexed sequence.
%
% For example:
%	For n=5, the indices of the aliased sequence are -2 -1 0 1 2.
%	Hence, lowIdx(n) = -2.
%	For n=4, the indices of the aliased sequence are -2 -1 0 1.
%	Hence, lowIdx(n) = -2.
%
% Yoel Shkolnisky 22/10/01

function idx = lowIdx(n)
	idx = -fix(n/2);
	