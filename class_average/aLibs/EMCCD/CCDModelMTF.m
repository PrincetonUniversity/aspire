function mtf=CCDModelMTF(f,iCamera)
% function mtf=CCDModelMTF(f)
% Return the analytical fit to the MTF of our CCD at 200kV
% as determined with MeasureMTF3.
% f is a frequency vector, with 0.5 (maximum value) being Nyquist.
% iCamera=1 : Yale F20 US4000
% iCamera=2 : Okazaki JEM 2200 TVIPS camera
% iCamera=3 : DE12
%
% f is a frequency vector of any dimension, with 0.5 (maximum value)
% being Nyquist.

if nargin<2
    iCamera=1;
end;

switch iCamera
    case 1
        f1=.0857;
        f2=.2250;
        f3=.3874;
        f=abs(f);
        mtf=(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));
        
    case 3
        f0=0.01513;
        a0=0.45748;
        f1=0.20599;
        f2=0.38772;
        f3=0.51004;
        mtf=a0./(1+(f/f0).^2)+(1-a0)*(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));

    case 4
        p=[0.026528     -0.21526     0.043531     0.099322      0.27167];
        f0=p(1); a0=p(2); f1=p(3); f2=p(4); f3=p(5);
        mtf=a0./(1+(f/f0).^2)+(1-a0)*(1+(f/f2).^2)./((1+(f/f1).^2).*(1+(f/f3).^3));
    
    otherwise
        error(['Invalid camera index ' num2str(iCamera)]);
end;

