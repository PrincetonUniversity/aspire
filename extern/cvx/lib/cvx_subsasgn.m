function y = cvx_subsasgn( y, varargin )
temp.type = '()';
temp.subs = varargin(1:end-1);
y = subsasgn( y, temp, varargin{end} );

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
