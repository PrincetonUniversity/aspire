function E=beadZ_params
% BEADZ_PARAMS  Parameters for a small Guassian on the z-axis.
% The position of the center of the Guassian corresponds to [0 0 1/2].
%
% See C1_params for the structure of E.
%
% Example:
%   vol=cryo_gaussian_phantom_3d('beadZ_params',65,1);
%   view3d(vol,0.5);
%
% Yoel Shkolnisy, November 2013.

E = [        
      1     1/12  1/12   1/12    0      0      1/2      0     0      0;      
     ];
