function [S cameraName]=CCDModelSpectrum(f,iCamera)
% Return the analytical fit to the spectrum of a CCD camera
% as determined with MeasureCCDSpectrum.  iCamera is an index into a table
% of fitted values.  At present the indices are
% iCamera=1 : Yale F20 US4000
% iCamera=2 : Okazaki JEM 2200 TVIPS camera
% iCamera=3 : DE12
%
% f is a frequency vector of any dimension, with 0.5 (maximum value)
% being Nyquist.

if nargin<2
    iCamera=1;
end;
% YaleF20 parameters...
%     g=[.1001
%         .2480
%         .3713
%         2.3783];
%     cameraName='YaleF20';


% Find our local directory
pa=fileparts(which('CCDModelSpectrum'));
% Retrieve parameters from SpectModel.mat in the local directory
if numel(pa)<1
    pa='.';
end;
load([pa '/SpectModel.mat']);

if iCamera>numel(models)
    error('CCDModelSpectrum: iCamera index too large');
end;
h=models(iCamera);

g=models(iCamera).spectPars;
aliasF0=models(iCamera).alias;
% aliasF0=0;%%%

cameraName=models(iCamera).camera;
if numel(g)<6
    S=g(4)./(1+(f/g(1)).^2+(f/g(2)).^4+(f/g(3)).^6);
else
    S=g(4)./(1+(f/g(1)).^2+(f/g(2)).^4+(f/g(3)).^6) + g(5)./(1+(f/g(6)).^2);
end;
if aliasF0>0
    af=aliasF0-f;
    S=S+g(4)./(1+(af/g(1)).^2+(af/g(2)).^4+(af/g(3)).^6) + g(5)./(1+(af/g(6)).^2);
end;

