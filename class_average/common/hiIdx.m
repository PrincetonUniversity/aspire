% function idx=hiIdx(n)
%
% Returns the maximal index for an aliased sequence of length n.
%
% n     the length of the indexed sequence.
%
% For example:
%   For n=5, the indices of the aliased sequence are -2 -1 0 1 2.
%   Hence, hiIdx(n) = 2.
%   For n=4, the indices of the aliased sequence are -2 -1 0 1.
%   Hence, hiIdx(n) = 1.
%
% Yoel Shkolnisky 22/10/01

function idx = hiIdx(n)
    idx = fix((n-0.5)/2);

% The above code actually performs:
%   if (mod(n,2)==0)
%       idx = n/2-1;
%   else    idx = (n-1)/2;
%   end
    
