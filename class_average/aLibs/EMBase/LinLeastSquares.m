function a=LinLeastSquares(F,y)
% function a=LinLeastSquares(F,y)
% Given a matrix F of function values (each column is values of 1 basis function)
% find the column vector of coefficients a that minimizes
% ||y-Fa||^2
% Thus each column of F corresponds to one of the functions to be fitted.
% The fitted function is obtain as F*a.

[n m]=size(F);
if numel(y)~=n
    error('Number of rows in F and y don''t match');
end;
y=y(:);

a=(F'*F)\(F'*y);
% a'

    