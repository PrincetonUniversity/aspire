function Q = randstiefel(varargin)
% RANDSTIEFEL Randomly draw a matrix from the Stiefel manifold
%
% Q = randstiefel(m,n) or randstiefel([m,n]) gives an m x min(m,n)
% orthogonal matrix drawn uniformly from the Stiefel manifold. 

[Q,R] = qr(randn(varargin{:}),0);

r = diag(R);
l = r./abs(r);

Q = Q*diag(l);
