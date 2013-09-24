function [ y, symm ] = cvx_s_banded( m, n, symm, lower, upper )

%CVX_S_BANDED (U,L)-banded matrices.

if nargin < 4,
    error( 'Bandwidth arguments missing.' );
end
if ~isnumeric( lower ) || length( lower ) ~= 1 || lower < 0 || lower ~= floor( lower ),
    error( 'Bandwidth arguments must be nonnegative integers.' );
elseif nargin < 5, 
    upper = lower;
elseif ~isnumeric( upper ) || length( upper ) ~= 1 || upper < 0 || upper ~= floor( upper ),
    error( 'Bandwidth arguments must be nonnegative integers.' );
end
if symm,
    lower = min( lower, upper );
    upper = 0;
end

c    = 0 : n - 1;
c    = c( ones( 1, m ), : );
r    = ( 0 : m - 1 )';
r    = r( :, ones( 1, n ) );
temp = r - c;
temp = temp <= lower & temp >= -upper;
r    = r( temp );
c    = c( temp );
nu   = length( r );
v    = ( 1 : nu )';

if symm,
    tt = r ~= c;
    r = [ r ; c(tt) ];
    c = [ c ; r(tt) ];
    v = [ v ; v(tt) ];
    symm = false;
end

y = sparse( v, r + m * c + 1, 1, nu, m * n );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
