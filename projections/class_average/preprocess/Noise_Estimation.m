function [ P ] = Noise_Estimation( projections )
%This function estimate noise power spectrum from the non-particle region
%   Input: projections
%   Output: P power spectrum
%   Zhizhen Zhao June 2013

%% Configure simulation parameters
verbose = 1;
norm_flag = 0;
biasflag = 1;
p = size(projections, 1); % sampling interval (image size)
K = size(projections, 3); % number of images
max_d = p*sqrt(2); % estimate length
a = 1; % correlation length
window_type = 'hann'; % check iwindow.m for details
radius_of_mask = floor(p/2)-1; % measurment radius > radius_of_mask
%% mask measurment window
I = cart2rad(p);
mask = zeros(p);
mask(I<radius_of_mask) = inf; % create a mask.
% mask(3:end-2, 3:end-2) = inf;
samples_idx = find(mask~=inf); % get points coordinate
%% measure 
fprintf('Measuring power spectrum from background pixels...\n');
[~, P] = epsdS_speedup( projections,samples_idx,biasflag, norm_flag, window_type,verbose );
P=real(P);
figure; imagesc(real(P));
end

