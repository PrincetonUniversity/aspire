function E=beadY_params
% BEADY_PARAMS  Parameters for a small Guassian on the y-axis.
% Note that the y axis is actually the first location parameter and not
% the second due to the flip of the x and y coordinates when using matrices
% as images in MATLAB.
% The position of the center of the Guassian corresponds to [0 1/2 0].
%
% See C1_params for the structure of E.
%
% Example:
%   vol=cryo_gaussian_phantom_3d('beadY_params',65,1);
%   view3d(vol,0.5);
%
% Yoel Shkolnisy, November 2013.

E = [        
      1     1/12  1/12   1/12    1/2      0     0      0     0      0;      
     ];
