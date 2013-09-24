function QuickLookSpectrum(image, pixA, name, T, ctf)
% function QuickLookSpectrum(image, pixA, name, T, ctf)
% Make a quick display of the image and its power spectrum.  Also shown is
% the radially-averaged spectrum.  The title is set to name, if given.
% The image argument can be a filename or a 2D array.  If it is an array, 
% the pixel size is given by pixA.  Otherwise the pixA value from the file,
% if available, is used.  The default pixA is 1.
% If no arguments are given, a file selector is put up.

exponent1=.3;
exponent2=.25;
neqbins=1e6;

fc=.3;
SetGrayscale;
if nargin < 5
    nrows=2;
else
    nrows=3;
end;
if nargin < 4
    T=1;
end;
if nargin < 3
    name='';
end;
if nargin < 2
    pixA=1;
end;

if nargin<1
    [image pa]=uigetfile('*','Select an image file');
    cd(pa);
end;

if ischar(image)  % must be a filename
    name=image;
    [m pixa ok]=ReadEMFile(image);
    if ~ok
        error('Invalid file type');
    end;
    if pixa~=0
        pixA=pixa;
    end;
    m=single(m);
%     m=RemoveOutliers(m);  % doesn't seem to affect the spectrum much.
else
    m=image;
end;

n=size(m);
ds=round(n(1)/512);
ndis=2*round(n/ds/2);
df=1/(pixA*n(1));

% Check if we need pre-whitening.
mc=Crop(m,ndis);
mc=mc-mean(mc(:));
smc=fftshift(abs(fftn(mc)).^2);
r=RadiusNorm(ndis);
mc1=(r>0.25 & r<=0.3);
mc2=(r>0.4 & r<=0.45);
s1=sum(mc1(:).*smc(:))/sum(mc1(:));
s2=sum(mc2(:).*smc(:))/sum(mc2(:));
PreWhiten = (s2/s1< 0.6 && s2/s1>.2);  % big roll-off
PreWhiten=0; %%%%%
if PreWhiten
     m=mePreWhiten(m);
     txt=' pre-whitened';
else
    txt='';
end;

subplot(nrows,2,1);
xrange=[0 pixA*n(1)]/10;
imacs(xrange,xrange,Downsample(m,ndis));
xlabel('nm');
title(name,'interpreter','none');
drawnow;

m=m-mean(m(:));
sp2=abs(fftn(m)).^2;
sp2(1,1)=0;
sp2d=abs(GaussFilt(BinImage(fftshift(sp2),ds),fc));
subplot(nrows,2,2);
frange=[-1 1]/(2*pixA);
% imacs(frange,frange,sp2d.^exponent);
imacs(frange,frange,ImEqualize((sp2d.^exponent1),neqbins)+0.3*imscale(sp2d.^exponent2));
xlabel('Spatial frequency, A^{-1}');
title([txt '  ' num2str(n) ' pixels']);
drawnow;

sp1=RadialPowerSpectrum(m);
nfreq=round(0.85*min(n)/2);
freqs=0:df:df*(nfreq-1);
hfAsymptote=mean(sp1(round(nfreq/2)+1:nfreq));
    Ts=hfAsymptote.*T.^2;
    if numel(Ts)<nfreq
        Ts(numel(Ts)+1:nfreq)=Ts(numel(Ts));
    end;
    Ts=Ts(:);
subplot(nrows,1,2);
fs=sqrt(freqs);
% save sp1 sp1  %%%
% sp2=sp1;
% load sp1
% sp2=load('sp1'); %%%%
semilogy((fs),[sp1(1:nfreq) Ts(1:nfreq)],'-');
% semilogy((fs),[sp1(1:nfreq) Ts(1:nfreq) sp2(1:nfreq)],'-');
% ,fs(1:11),sp1(1:11),'.');
SquareRootTicks(fs);
xlabel('Spatial frequency, A^{-1}');
title(s2/s1);

if nargin>4  % plot the ctf too.
subplot(nrows,1,3);
plot(fs,ctf(1:nfreq));
SquareRootTicks(fs);
xlabel('Spatial frequency, A^{-1}');
end;
drawnow;
