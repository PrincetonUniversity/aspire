function ShowSpectrum(m,res,fraction, exponent)
% function ShowSpectrum(m,res,exponent,fraction)
% Show the central fraction of the power spectrum of m,
% where the grayscale represents the spectral values raised to the
% exponent.  res is the scale in A/pixel.
% Default values are fraction=1/4, exponent=1/4.

if nargin<4
    exponent=0.25;
end;
if nargin<3
    fraction=.25;
end;
n=size(m,1);
m=m-mean(m(:));
am=fftshift(abs(fftn(m)));
n1=2*round(fraction*n/2);
am=Crop(am,n1);

SetGrayscale;
imacs(am.^(2*exponent));

df=10/(res*n);  % Frequency increment per pixel
LabelImage(n1,df);
xlabel('Spatial frequency, nm^{-1}');
