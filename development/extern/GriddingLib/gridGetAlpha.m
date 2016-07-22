function alpha=gridAlpha(kernelsize);
% Return an alpha value for gridding interpolation.  Invalid arguments
% cause 0 to be returned.
%
    alphavals=[0 0 5.2 0 10.2 0 13 0 17];  % nice alpha values for the 1.25 x oversampling.
if kernelsize>9
    alpha=0;
else
    alpha=alphavals(kernelsize);
end;
if alpha==0
    error(['invalid kernelsize ' num2str(kernelsize)]);
end;
