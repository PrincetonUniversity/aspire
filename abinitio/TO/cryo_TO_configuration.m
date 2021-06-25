function [n_theta, n_r, rmax, resolution, viewing_angle, inplane_rot_degree] = cryo_TO_configuration()

% configuration file.
%       
%   Written by Adi Shasha January 2021. 

% image parameters
n_theta = 360;                                      % L radial lines.
n_r     = 89;                                       % Number of equispaced samples along each radial line.
rmax    = 2;  

% generating candidates set parameters
resolution         = 75;                            % the number of samples per 2*pi (see genRotationsGrid for more details).
viewing_angle      = 0.996;                         % the viewing angle threshold.
inplane_rot_degree = 5;                             % the inplane rotation degree threshold.