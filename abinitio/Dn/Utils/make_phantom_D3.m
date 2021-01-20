function volD3=make_phantom_D3
%
% Generate a phantom with D3 symmetry.
% The phantom is accurate to 10^-15.
%
% Example:
%   view3d(volD3);
%
% Yoel Shkolnisky, November 2020.

% Generate an off center Gaussian
vol=cryo_gaussian_phantom_3d('football_params',89,1);

% Apply Cn group elements
vol120=fastrotate3z(vol,120);
vol240=fastrotate3z(vol,240);
volz=vol+vol120+vol240; % At this point the volume volz is C3.

% Apply 180 degrees rotation
volx=fastrotate3x(volz,180);

% Combine both volumes
volD3=volx+volz;

% Test code. All errors should be small
err=volD3-fastrotate3z(volD3,120);
assert(norm(err(:)/norm(volD3(:)))<1.0e-14)
err=volD3-fastrotate3z(volD3,240);
assert(norm(err(:)/norm(volD3(:)))<1.0e-14)
err=volD3-fastrotate3x(volD3,180);
assert(norm(err(:)/norm(volD3(:)))<1.0e-14)

