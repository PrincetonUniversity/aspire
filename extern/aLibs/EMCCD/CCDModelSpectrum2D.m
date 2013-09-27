function [S cameraName]=CCDModelSpectrum2D(iCamera,n)
% Return the analytical fit to the spectrum of a CCD camera
% as determined with MeasureCCDSpectrum2D.  iCamera is an index into a table
% of fitted values.  At present the indices are
% iCamera=1 : Yale F20 US4000
% iCamera=2 : Okazaki JEM 2200 TVIPS camera
% iCamera=3 : DE12
% iCamera=4 : BrandeisFalcon1
% (For all but iCamera=3, we just call CCDModelSpectrum with a 2D array of
% frequencies)
% n is a scalar or 2-element vector giving the size of the image to be
% returned. Default is the full camera size (3k x 4k for DE12)
% f=0 is in the center, and the edges represent f=0.5 (Nyquist)

if nargin<2
    switch iCamera
        case 3
            n=[3072 4096];
        otherwise
            n=[4096 4096]; % indices 1,2,4
    end;
end;
if numel(n)==1 % square image
    n=[n n];
end;
% YaleF20 parameters...
%     g=[.1001
%         .2480
%         .3713
%         2.3783];
%     cameraName='YaleF20';


% Find our local directory
pa=fileparts(which('CCDModelSpectrum2D'));
% Retrieve parameters from SpectModel2D.mat in the local directory
if numel(pa)<1
    pa='.';
end;
load([pa '/SpectModel2D.mat']);

if iCamera>numel(models)
    error('CCDModelSpectrum2D: iCamera index too large');
end;

p=models(iCamera).spectPars;
cameraName=models(iCamera).camera;

if iCamera==3   % special case of the DE camera
    
    [f theta]=RadiusNorm(n);
    % Aliasing code copied from ccFit2DSpectrum
    ang=max(abs(cos(theta)),abs(sin(theta)));
    fAliased=1./ang-f;
    y0=ccDE12Model(f,theta,p);
    ya=ccDE12Model(fAliased,theta,p);
    S=y0+ya;
    
else
    f=RadiusNorm(n);
    S=CCDModelSpectrum(f,iCamera);
end;
