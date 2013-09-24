function cvx_optval = prod_inv( x, dim, p )

%DET_INV   Internal cvx version.

error( nargchk( 1, 3, nargin ) ); %#ok
if ~isreal( x ), 
    error( 'First argument must be real.' ); 
end
sx = size( x );
if nargin < 2 || isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim ),
    error( 'Second argument must be a positive integer.' );
end
sx( end + 1 : dim ) = 1;
if nargin < 2,
    p = 1;
elseif ~isnumeric( p ) || ~isreal( p ) || numel( p ) ~=  1 || p <= 0,
    error( 'Third argument must be a positive scalar.' );
end

if cvx_isconstant( x ),
    
    cvx_optval = cvx( prod_inv( cvx_constant( x ), dim, p ) );

elseif sx( dim ) == 1,
    
    cvx_optval = inv_pos( x );
    
else

    sy = sx;
    sy( dim ) = 1;
    y = [];
    cvx_begin
        epigraph variable y( sy )
        geo_mean( cat( dim, x, y ), dim, [ ones(n,1) ; p ] ) >= 1; %#ok
    cvx_end

end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
