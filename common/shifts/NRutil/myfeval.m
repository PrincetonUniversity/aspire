function y=myfeval(f,x,data)
%
% Evaluate the function f at the point x using the additional data.
% If data is an empty array, then the function is evaluated using x only.
%
% Yoel Shkolnisky, July 2008.

if isempty(data)
    y=feval(f,x);
else
    y=feval(f,x,data{:});
end
