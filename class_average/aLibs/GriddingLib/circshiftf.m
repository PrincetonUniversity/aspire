function ms=circshiftf(m,shiftVector)
% function ms=circshiftf(m,shiftVector)
% 2D shift by multiplying the FT by the appropriate complex exponential.

% % Example test code
% subplot(2,2,1);
% m=fuzzymask(128,2,10,0.2);
% shiftVector=[10.5 0];
% ms=FourierShift(m,shiftVector);
% subplot(2,2,2);
% imacs(ms+m);
% drawnow;


mf=fftn(m);
[n n1]=size(m);
[X Y]=ndgrid(-n/2:n/2-1);
V=shiftVector(1)*X+shiftVector(2)*Y;
mfs=mf.*exp(-1j*fftshift(V)*2*pi/n);
ms=real(ifftn(mfs));

