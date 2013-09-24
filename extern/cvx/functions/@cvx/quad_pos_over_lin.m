function cvx_optval = quad_pos_over_lin( x, y, varargin ) %#ok

%QUAD_POS_OVER_LIN   Internal cvx version.

error( nargchk( 2, 3, nargin ) ); %#ok
if ~isreal( x ), 
    error( 'First input must be real.' ); 
end
x2 = [];
cvx_begin
    variable x2( size(x) )
    minimize quad_over_lin( x2, y, varargin{:} );
    x2 >= x; %#ok
cvx_end

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

