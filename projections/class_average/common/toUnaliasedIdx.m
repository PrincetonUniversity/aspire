% function y=toUnaliasedIdx(idx,n)
%
% Converts an index from the range -n/2...n/2-1 to an index in the range 1...n.
% Both odd and even values of n are handled. 
%
% idx    An index from the range -n/2...n/2-1.
% n      The range of indices.
%
% Returns the index "idx" scaled to the range 1...n.
% 
% For example:
%       n   = 5;
%       idx = -1;
%       toUnaliasedIdx(idx,n) will return 2:
%           -2 -1 0 1 2
%            ^
%        the index of -1 is 2 if scaled to 1...n.
%
% Yoel Shkolnisky 22/10/01

function y=toUnaliasedIdx(idx,n)

% verify that idx and n are scalars
if (prod(size(idx))~=1) | (prod(size(n))~=1)
    error('idx and n must be scalars');
end

y=idx+floor(n/2)+1;

%Revision Record
%  8/1/03 Yoel Shkolnisky    Optimization
%     previous code:    
%          if (mod(n,2)==0)
%              y=idx+n/2+1;
%           else y=idx+(n-1)/2+1;
%          end
%     Using floor (built-in function) is much faster.
