function [ P,R,R2 ] = cryo_noise_estimation( projections,radius_of_mask )
%
% CRYO_NOISE_ESTIMATION Estimate 2D noise power spectrum of the projections.
%
% [ P,R, R2 ] = cryo_noise_estimation( projections )
%   Estimate the 2D isotropic power spectrum of the noise in the stack
%   "projections" using in each image only the pixels at the corners.
%
% Input parameters:
%   projections     Stack of projections. imstack(:,:,k) is the k'th
%       projection in the stack. Projections must be square and can have
%       odd or even dimensions.
%
% Output parameters:
%   P   2D power spectrum function. If each image is of size pxp, then P2
%   is of size (2p-1)x(2p-1). P is always real.
%   R   1D isotropic autocorrelation function.
%   R2  2D isotropic autocorrelation function.

% Zhizhen Zhao June 2013
% Revised: Yoel Shkolnisky, October 2014, July 2015.
% Y.S. March 3, 2016    Add R2 as an output variable.


verbose = 1;
p = size(projections, 1); % sampling interval (image size)

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
[P,R, R2] = cryo_epsdS( projections,noise_idx,floor(p/3),verbose );
P=real(P);
end
