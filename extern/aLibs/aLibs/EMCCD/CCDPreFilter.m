function h=CCDPreFilter(n,iCamera,fc,binning)
% function h=CCDPreFilter(n,iCamera,fc,binning)
% Compute pre-whitening filter for a CCD camera.  The returned filter
% function has zero frequency in the center, so it must be fftshifted.
% The optional argument
% iCamera is the index of the camera (1=YaleF20; 2=NIPS2200).  Default is 1.
% If fc is given, it is the corner frequency relative to (Nyquist/binning)
% of an erf filter; the default value is 0.95.
% The argument binning is also optional; default is 1.  If the image has
% been binned--decreasing the maximum frequency to (f_Nyquist/binning)--
% then give this value.
%
%  Here is an example of prewhitening image m to make the filtered image mf:
%  n=size(m,1);
%  h=CCDPreFilter(n,1);  % Yale F20 camera selected
%  mf=real(ifftn(fftn(m).*fftshift(h));

if nargin<2
    iCamera=1;
end;
if nargin<3
    fc=.975;
end;
if nargin<4
    binning=1;
end;
f=Radius(n)/(n*binning);

h=1./sqrt(CCDModelSpectrum(f,iCamera));
if fc>0
    h=h.*fuzzymask(n,2,fc*n/2,n*(1-fc));
end;
