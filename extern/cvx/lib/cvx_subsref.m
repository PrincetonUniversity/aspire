function y = cvx_subsref( x, varargin )
temp.type = '()';
temp.subs = varargin;
y = subsref( x, temp );

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
