function [ P ] = Noise_Estimation( projections )
%This function estimate noise power spectrum from the non-particle region
%   Input: projections
%   Output: P power spectrum
%   Zhizhen Zhao June 2013
% Revised: Yoel Shkolnisky, October 2014.

% addpath /u/zhizhenz/matlab_cryo/epsdS
%% Configure simulation parameters
verbose = 1;
p = size(projections, 1); % sampling interval (image size)
radius_of_mask = floor(p/2)-1; % measurment radius > radius_of_mask
%% mask measurment window
I = cart2rad(p);
mask = zeros(p);
mask(I<radius_of_mask) = -1; % create a mask.
samples_idx = find(mask~=-1); % get points coordinate
%% measure 
fprintf('Measuring power spectrum from background pixels...\n');
P = cryo_epsdS( projections,samples_idx,floor(p/3),verbose );
P=real(P);
end

