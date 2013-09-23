function [b,m,q]=nufftparams(precision)
%
% Set the parameters for the non-equally spaced FFT according to the
% required precision.
%
% Yoel Shkolnisky, December 2009.

if (~strcmpi(precision,'single')) && (~strcmpi(precision,'double'))
    precision='single';
    warning('MATLAB:unkownOption','Unrecognized precsion. Using ''single''.');
end

if strcmpi(precision,'double')
    b=1.5629;
    m=2;
    q=28;
else %single precision
    b=0.5993;
    m=2;
    q=10;
end

if mod(q,2)==1 % Just for safety.
    error('Set q to an even integer');
end
