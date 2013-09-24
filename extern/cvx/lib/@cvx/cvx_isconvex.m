function y = cvx_isconvex( x, full )
error( nargchk( 1, 2, nargin ) );
y = cvx_vexity( x );
if nargin < 2,
    y = nnz( y < 0 ) == 0;
else
    y = cvx_reshape( y >= 0, x.size_ );
end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
