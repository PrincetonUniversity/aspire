function [w, ov]=gridMakeKaiserTable(kernelsize,mode)
% Create the interpolated Kaiser or sinc-Kaiser kernel, with oversampling factor ov
% The mode variable is either 'grid' or 'sinc'.
% the returned w array is ov x kernelsize in size.

persistent w1;  % Cache for the table

nw=kernelsize;
mo=lower(mode(1));

switch lower(mode(1))
    case 'g' % gridding
        sincmode=0;
    case 's' % sinc
        sincmode=1;
    otherwise
        error(['invalid mode:' mode]);
end;

if mo=='s' % sinc mode
    sincmode=1;
    % parameters and static variables
    ov=128;  % oversampling of window for look-up
    % alphavals=[0 0 2.5 0 4.5 0 5.7 0 6.8]; % values for minimum aliasing.
    alphavals=[0 0 3.5 0 5 0 6.8 0 6.8]; % these values seem to work better.
    alpha=alphavals(nw);
    if alpha==0
        error('Invalid kernel size; it must be odd.');
    end;

elseif mo=='g' % gridding: use kaiser kernel
    sincmode=0;
    alpha=gridGetAlpha(kernelsize);
    ov=1024;  % oversampling of window for look-up
else
    error(['invalid mode:' mode]);
end;

if (numel(w1)~=nw*ov) % if the cached table is absent or not the right size
    % Create the interpolation kernel
    epsi=1e-8;  % Crude way to avoid divide by zero in the sinc function.
    dw=1/ov;  % fractional increment.  We assume ov is even
    k=(-nw/2+dw/2:dw:nw/2-dw/2)';  % nw*ov space points
    w=kaiser(nw/2,alpha,k);  % Compute the Kaiser window
    if sincmode
        w=w.*sin(pi*k+epsi)./(pi*k+epsi);  % multiply by sinc function
    end;
    w=w*ov/sum(w);  % normalize it.
    w1=zeros(ov,nw);
    % Make the 1D lookup table w1(i,:).  i=1 corresponds to a shift between -0.5 and
    % 0.5+dw, so we assign it the mean shift -0.5+dw/2.
    % i=ov corresponds to a mean shift of 0.5-dw/2.
    for i=1:ov
        w1(i,:)=w(ov-i+1:ov:nw*ov);
    end;
end;
w=w1'; % transpose to (nw x ov) size.
