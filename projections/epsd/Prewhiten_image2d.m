function [ proj ] = Prewhiten_image2d(proj, noise_response )
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, May 2016

warning('This function is depreacted. Call cryo_prewhiten instead.');
[ proj ] = cryo_prewhiten(proj, noise_response );