function a = in( x, y )

error( nargchk( 2, 2, nargin ) );
b = newcnstr( evalin( 'caller', 'cvx_problem', '[]' ), x, y, '==' );
if nargout, a = b; end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
