function y = cvx_isconcave( x, full ) %#ok
error( nargchk( 1, 2, nargin ) ); %#ok
if nargin == 2,
    y = ~imag( x );
elseif isreal( x ),
    y = true;
else
    y = nnz(imag(x)) == 0;
end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
