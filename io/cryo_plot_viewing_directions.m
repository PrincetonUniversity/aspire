function [h1,h2]=cryo_plot_viewing_directions(rotations)
%
% CRYO_PLOT_VIEWING_DIRECTIONS  Plot the distribution of the viewing 
%                               directions. 
%
% 
% cryo_plot_viewing_directions(rotations)
%   "rotations" are the estimated rotations as returned, for example, by
%   the synchronization algorithm. Returns a handle to the generated
%   figure.
%
% Yoel Shkolnisky, March 2015.

viewingdirs=rotations(:,3,:); % Viewing direction is the third column of 
                              % each rotation matrix.
viewingdirs=squeeze(viewingdirs);
thetas=atan2d(viewingdirs(2,:),viewingdirs(1,:));
%phis=acosd(viewingdirs(3,:));
rhos=sqrt(1-viewingdirs(3,:).^2);

h1=figure;
%scatter(thetas(:),phis(:));
polar(thetas/180*pi,rhos,'.');

h2=figure;
scatter3(viewingdirs(1,:),viewingdirs(2,:),viewingdirs(3,:));