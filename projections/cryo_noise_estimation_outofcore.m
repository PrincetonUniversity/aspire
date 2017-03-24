function [ P,R,R2 ] = cryo_noise_estimation_outofcore( instackname,radius_of_mask )
%
% CRYO_NOISE_ESTIMATION_OUTOFCORE Estimate 2D noise power spectrum of the projections.
%
% [ P,R, R2 ] = cryo_noise_estimation( instackname, radius_of_mask )
%   Estimate the 2D isotropic power spectrum of the noise in the
%   projections of the MRC file named instackname, using in each image only
%   the pixels whose radius is large than radius_of_mask.
%   The function does not load the entire stack into memory.
%   
% See cryo_noise_estimation for more details.
%
% Example:
%   cryo_noise_estimation_outofcore('instack.mrc');
%   cryo_noise_estimation_outofcore('instack.mrc',32);
%
% Yoel Shkolnisky, May 2016.

verbose = 1;
instack=imagestackReader(instackname);
p = instack.dim(1); % sampling interval (image size)

if nargin<2    
    radius_of_mask = floor(p/2)-1; % measurment radius > radius_of_mask
end

%% Mask measurment window
I = cart2rad(p);
mask = zeros(p);
mask(I<radius_of_mask) = -1;  % Create a mask.
noise_idx = find(mask~=-1);   % Get points coordinate

%% Estimate power spectrum.
%log_message('Measuring power spectrum from background pixels...');
[P,R, R2] = cryo_epsdS_outofcore( instackname,noise_idx,floor(p/3),verbose );
P=real(P);
end

