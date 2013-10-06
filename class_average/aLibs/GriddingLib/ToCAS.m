function r=ToCAS(c)
% function r=ToCAS(c)
% Convert the complex vector or matrix c to a corresponding real one.
% The dimension(s) of c must be even
r=real(c)+imag(c);
