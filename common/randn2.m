function [ output ] = randn2(varargin)
% randn2  Calls rand and applies inverse transform sampling to the output.
 
output = rand(varargin{:});
output = sqrt(2) * erfinv(2 * output - 1);
 
end