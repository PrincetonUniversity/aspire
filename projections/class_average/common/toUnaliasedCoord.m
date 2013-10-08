%function uCoord = toUnaliasedCoord(aCoord,N)
%
% Converts indices from the range -n/2...n/2-1 to indices in the range 1...n.
% Both odd and even values of n are handled. 
%
% The functions accepts a vector of aliased indices (aCoord) and a vector of ranges
% (N) from which the indices are taken. It converts each aliased index aCoord(k) into an
% unaliased index uCoord(k) in the range 1...N(k). If the vector of ranges
% (N) is a scalar, then the function assumes that all indices are taken
% from the same range N. This allows calling the function on a vector of
% indices:
%       toUnaliasedCoord([-1 1 2],5)
% instead of
%       toUnaliasedCoord([-1 1 2],[5 5 5])
%
% Input:
%   aCoord    Vector of aliased indices. Must be 1-D row vector.
%   N         Vector that contains the range of each index. Must be 1-D row
%             vector or a scalar. If N is scalar it is used for all
%             coordinates.
% Output:
%   uCoord    Vector of unaliased indices.
%
% If N is not a scaler, the vectors aCoord and N must have the same length.
%
% Yoel Shkolnisky 8/1/03

function uCoord = toUnaliasedCoord(aCoord,N)

% Remove validity test in increase performance
if (length(size(aCoord))~=2) | (length(size(N))~=2)
   error('Input parameters must be 1D row vectors');
end

if (size(aCoord,1)~=1) | (size(N,1)~=1)
   error('Input parameters must be 1D row vectors');
end

if (length(N)~=1) & (length(aCoord)~=length(N))
   error('N must be scalar or length of coordinates vector and ranges vector must be the same');
end

% The boolean flag scalar is 1 if the ranges vector N is a scalar.
scalar=0;
if length(N)==1
    scalar=1;
end


uCoord = cell(size(aCoord));
for k=1:length(aCoord)
    if scalar
        uCoord{k} = toUnaliasedIdx(aCoord(k),N);
    else
        uCoord{k} = toUnaliasedIdx(aCoord(k),N(k));
    end
end
