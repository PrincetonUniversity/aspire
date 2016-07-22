function imacx(m,power)
%  imacx(m), imacx(m,power)
% Autoscaled plot of complex-valued matrix m.  Intensity of the plot is
% determined by abs(m), and the color represents the phase.
% If the argument 'power' is given, then the intensity is abs(m).^power.
% Typically power<1 to expand the dynamic range of the display.

m=squeeze(m);

if nargin<2
    r=abs(m);
else
    r=abs(m).^power;
end;

t=atan2(imag(m),real(m));  % theta
mx=max(max(r));
if mx > 0
    r=r/mx;
end;

% wrap t about zero, not -pi
td=t/(2*pi)+(t<0);
d=64*(1-eps); % magic factor
md=257+floor(r*d)+64*floor(d*td);

% draw the image
image(md');
axis xy
