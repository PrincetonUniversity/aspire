% function g=slowChirp(x,A,W,M)
%
% Chirp Z-transform of the sequence x on the contour defined by
% A*W^(-k), where k is a sequence of M indices centered around the origin.
% For example, for M=5, the sequence k is -2,-1,0,1,2. For M=4, the sequence 
% k is -2,-1,0,1.
% The chirp Z-transform is computed directly using O(n^2) operations.
%
% x    The sequence whose chirp Z-transform should be computed. Can be of odd or even length.
% A    Arbitrary complex number.
% W    The complex ratio between two consecutive points on the contour.
% M    The length of the output sequence. If not specified, the default value 
%      is M=length(x);
%
% Returns the chirp Z-transform of the sequence x define by
%                n/2-1
%          X(Z) = sum  x(j)Z^(-j)
%                j=-n/2
% along the contour Z_k = AW^(-k)     k=-M/2...M/2-1.
%
% For example, for x = [1 2 3 4 5], the call
%     slowChirp(x,1,exp(-2*pi*i/5),5)
% computes the aliased DFT of x.
% 
% Yoel Shkolnisky 6/1/03

function g=slowChirp(x,A,W,M)

if nargin<4
   M=length(x);
end

n=length(x);
g=zeros(1,M);

for k=lowIdx(M):hiIdx(M)
   Z = A*W^(-k);
   acc = 0;
   for j=lowIdx(n):hiIdx(n)
      acc = acc + x(toUnaliasedIdx(j,n))* Z^(-j);
   end;
   g(toUnaliasedIdx(k,M)) = acc;
end
