function [ P,R ] = Noise_Estimation( projections )
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, August 2015

warning('This function is depreacted. Call cryo_noise_estimation instead.');
[ P,R ] = cryo_noise_estimation( projections );