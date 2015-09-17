function h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
                SphericalAberration,pixA,AmplitudeContrast)
% CTF   Compute constract transfer fuuction
%
% h=CTF(n,voltage,DefocusU,DefocusV,DefocusAngle,...
%                SphericalAberration,pixA,AmplitudeContrast)
%
% Compute the contrast transfer function corresponding an n x n image with
% the sampling interval DetectorPixelSize.
%
% Input parameters:
%  n                   Size of the image (nxn). 
%  voltage             Microscope electric potential (in Kvolts).
%  DefocusU,DefocusV   Two defocus values describing the defocus in two
%                      perpendicular directions in an image when
%                      astigmatism is present (in nm).
%  DefocusAngle        Angle between the first direction (described by DF1)
%                      and the X-axis (in radians).
%  SphericalAberration Spherical aberration, typically denoted by Cs 
%                      (in mm).
%  pixA                Pixel size (in Angstrom). This is computed by
%                      dividing the detector pixel size (something like 14
%                      microns) by the magnification.
%  AmplitudeContrast   Amplitdue contrast.
%
% Output:
%  h    Array of size nxn. Zero frequnecy amplitdue is given at
%       h(n/2+1,n/2+1). Thus you must use ifftshift() on the result before
%       performing a Fourier transform.  For example, to simulate the
%       effect of the CTF on an image m, do this: 
%           fm=fftshift(fft2(m));
%           cm=ifft2(ifftshift(fm.*ctf()));
% 
% Algorithm implemented according to
% Mindell, J. A.; Grigorieff, N. (2003). "Accurate determination of local
% defocus and specimen tilt in electron microscopy". Journal of structural
% biology 142 (3): 334�347. PMID 12781660  
%
% For more information see
% http://i2pc.cnb.csic.es/emx/LoadDictionaryFormat.htm?type=Convention 
%
% Example:
%   n=256;
%   V=300;
%   DF1=2.334469921900000e+03;
%   DF2=2.344594921900000e+03;
%   theta= 36.700001;
%   Cs=2.0;
%   pixA=1.4;
%   A=0.1;
%   h=cryo_CTF_Relion(n,V,DF1,DF2,theta,Cs,pixA,A);
%   imagesc(h);
%
% Yoel Shkolnisky, July 2014.

lambda=1.22639./sqrt(voltage*1000+0.97845*voltage.^2);  % Wavelength in nm.

BW=1/(pixA/10); % Divide by 10 to make pixel size in nm. BW is the 
                % bandwidth of the signal corresponding to the given pixel
                % size. 
                             
[s, theta]=RadiusNorm(n,fctr(n));
s=s.*BW; % RadiusNorm returns radii such that when multiplied by the 
         % bandwidth of the signal, we get the correct radial frequnecies
         % corresponding to each pixel in our nxn grid. 

DFavg=(DefocusU+DefocusV)/2;
DFdiff=(DefocusU-DefocusV);
df=DFavg+DFdiff*cos(2*(theta-DefocusAngle))/2;

k2=pi.*lambda.*df;
k4=pi/2*10^6*SphericalAberration*lambda^3; %10^6 converts Cs from mm to nm.
chi=k4.*s.^4-k2.*s.^2;
h=sqrt(1-AmplitudeContrast^2).*sin(chi)-AmplitudeContrast.*cos(chi);

% h=sqrt(1-AmplitudeContrast^2).*sin(k4.*s.^4-k2.*s.^2)-...
%     AmplitudeContrast.*cos(k4.*s.^4-k2.*s.^2);

