% function g=OptimizedChirpZ(x,A,W,M)
%
% Optimized implementation of the function chirpZ.
%
% The function saves precomputed factors between consecutive executions.
% Therefore, if the function is called many times for different sequences x
% with the same parameters A,W,M, the function will execute faster than the
% function chirpZ.
%
% See chirpZ for more information
%
% Yoel Shkolnisky 6/1/03 

function g=OptimizedChirpZ(x,A,W,M)

% LAST_CHIRP_Z_EXECUTION_RESULTS is a persistent cache that contain all data from
% the previous execution that can be pre-computed.
persistent LAST_CHIRP_Z_EXECUTION_RESULTS;

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

% check if we can use the cache
if isempty(LAST_CHIRP_Z_EXECUTION_RESULTS)
   LAST_CHIRP_Z_EXECUTION_RESULTS.n = 0;
   LAST_CHIRP_Z_EXECUTION_RESULTS.M = 0;
   LAST_CHIRP_Z_EXECUTION_RESULTS.W = 0;
end

if (LAST_CHIRP_Z_EXECUTION_RESULTS.n ~= n) | (LAST_CHIRP_Z_EXECUTION_RESULTS.M ~= M) | (LAST_CHIRP_Z_EXECUTION_RESULTS.W ~= W)
   y_factors = W.^((j.^2)/2);
   v=W.^(-(j2.^2)/2);
   V=cfft(v);
   outIdx = lowIdx(M):hiIdx(M);     
   post_factors = W.^(outIdx.^2/2);
   n_range = toUnaliasedIdx(lowIdx(n),pl):toUnaliasedIdx(hiIdx(n),pl);
   M_range = toUnaliasedIdx(lowIdx(M),pl):toUnaliasedIdx(hiIdx(M),pl);
   %update the cache
   LAST_CHIRP_Z_EXECUTION_RESULTS = struct('n',n,'M',M,'W',W,'y_factors',y_factors,'V',V,'post_factors',post_factors,'n_range',n_range,'M_range',M_range);
%   fprintf('MISS\n');
%else
%   fprintf('HIT\n');
end   

y(LAST_CHIRP_Z_EXECUTION_RESULTS.n_range) =x.*(A.^(-j)).*(LAST_CHIRP_Z_EXECUTION_RESULTS.y_factors);
% Convolve the arrays y and v
Y=cfft(y);
G=Y.*LAST_CHIRP_Z_EXECUTION_RESULTS.V;
g=icfft(G);
% Extract relevant portion of the array - the portion the corresponds to -n/2<=k<=n/2
g=g(LAST_CHIRP_Z_EXECUTION_RESULTS.M_range);
% Postmultiplication
g=g.*(LAST_CHIRP_Z_EXECUTION_RESULTS.post_factors);
