function cvx_optval = square_pos( x )

%SQUARE_POS   Internal cvx version.

error( nargchk( 1, 1, nargin ) ); %#ok
cvx_optval = pow_pos( x, 2 );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
