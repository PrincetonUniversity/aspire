function [ proj, filter, nzidx] = Prewhiten_image2d(proj, noise_response )
% Pre-whiten a stack of projections using the power spectrum of the noise.
%  
%  Input:
%   proj            Stack of images.
%   noise_response  2d image with the power spectrum of the noise.
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, May 2016

warning('This function is depreacted. Call cryo_prewhiten instead.');
[ proj, filter, nzidx] = cryo_prewhiten(proj, noise_response );
