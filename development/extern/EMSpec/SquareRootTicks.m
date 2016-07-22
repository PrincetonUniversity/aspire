function SquareRootTicks(xvals)
% function SquareRootTicks(xvals,reciprocal)
% Put x-axis ticks in the appropriate place for plotting square-root
% values.
% e.g. plot(sqrt(freqs),spect); SquareRootTicks(freqs);
% if nargin<2
% reciprocal=0;
% end;
xmin=min(xvals);
xmax=max(xvals);
dx=Step125(xmax/20);
dx2=Step125(dx/2);
dx4=Step125(dx2/2);
x0=floor(xmin/dx)*dx;
tickVals=[x0:dx4:dx2-dx4 dx2:dx2:dx-dx2 dx:dx:xmax+dx-eps]';
% tickVals
% if reciprocal
%     vals=1/tickVals;
% else
    vals=tickVals;
% end;
labels=num2str(vals);
set(gca,'xtick',sqrt(tickVals),'xticklabel',labels);
