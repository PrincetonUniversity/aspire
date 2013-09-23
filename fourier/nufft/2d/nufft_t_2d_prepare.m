function precomp=nufft_t_2d_prepare(x,n,precision)
%
% Prepare interpolation tables for use with nufft_t_2d_execute.
%
% Input parameters
%    x   Sampling points in Fourier space. An array with two columns of
%        real numbers in the range [-pi,pi].  
%    n   Size of the transformed image (as will be provided later to the
%        function nufft_t_v3(beta,precomp)
%    precision  'double' or 'single'. Default is 'single'.
%
% Output
%    precomp   Data structure to be passed to nufft_t_v3(beta,precomp)
%              whenever the the size of the transformed image and the
%              sampling points in the frequency domain have not been
%              changed.
% 
% Example:
%    precomp=nufft_t_prepare_v3(freqs,64,'single');
%    pf=nufft_t_v3(im,precomp);
%
% Yoel Shkolnisky, January 2008
%
% Revisions:
% Filename changed from nufft_t_prepare_v3.m to nufft_t_2d_prepare. Minor
% code fixes. (Y.S. December 22, 2009).

if nargin<3
    precision='single';
end

if (~strcmpi(precision,'single')) && (~strcmpi(precision,'double'))
    precision='single';
    warning('GCAR:malformedArgument','Unrecognized precsion. Using ''single''.');
end

if strcmpi(precision,'double')
    b=1.5629;
    m=2;
    q=28;
else %single precision
    b=0.5993;
    m=2;
    q=10;
end

if length(size(x))~=2
    error('x must be a mx2 array')
end

if size(x,2)~=2
    error('x must be a mx2 array');
end

nu=round(x*m*n/(2*pi));
precomp=nufftt2dpreparemx(x,nu,n,m,b,q);
   
