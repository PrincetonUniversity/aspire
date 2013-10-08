% function g=ChirpZ(x,A,W,M)
%
% Chrip Z-transform of the sequence x on the contour defined by
% A*W^(-k), where k is a sequence of M indices centered around the origin.
% For example, for M=5, the sequence k is -2,-1,0,1,2. For M=4, the sequence 
% k is -2,-1,0,1.
%
% The chirp Z-transform is computed using O(nlogn) operations.
%
% x    The sequence whose chirp Z-transform should be computed. Can be of odd or even length.
% A    Arbitrary complex number.
% W    The complex ratio between two consecutive points on the contour.
% M    The length of the output sequence. If not specified, the default value
%      is M=length(x);
%
% Returns the chirp Z-transform of the sequence x define by
%                n/2-1
%          X(Z) = sum x(j)Z^(-j)
%                j=-n/2
% along the contour Z_k = AW^(-k)     k=-M/2...M/2-1.
%
% Yoel Shkolnisky 6/1/03

function g=ChirpZ(x,A,W,M)

if nargin<4
   M=length(x);
end

n=length(x);
j=lowIdx(n):hiIdx(n);     
pl = M+n; % the required total length for the convolution.
j2= lowIdx(pl):hiIdx(pl);

% Create the array y of length pl and place the terms of the sequence, defined 
% in the paper, in the middle (centered about zero).
x = x(:).'; % ensure that x is a row vector
y=zeros(1,pl);
y(toUnaliasedIdx(lowIdx(n),pl):toUnaliasedIdx(hiIdx(n),pl)) =x.*(A.^(-j)).*(W.^((j.^2)/2));

% Create the array v
v=W.^(-(j2.^2)/2);

% Convolve the arrays y and v
Y=cfft(y);
V=cfft(v);
G=Y.*V;
g=icfft(G);

% Extract relevant portion of the array - the portion the corresponds to -n/2<=k<=n/2
g=g(toUnaliasedIdx(lowIdx(M),pl):toUnaliasedIdx(hiIdx(M),pl));

% Postmultiplication
outIdx = lowIdx(M):hiIdx(M);     
g=g.*(W.^(outIdx.^2/2));
