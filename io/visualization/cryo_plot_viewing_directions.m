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
%   Each point correspond to the (theta,phi) coordintaes of the view
%   direction, colored by the polar angle phi (restricted to [0,90]).
%
% Yoel Shkolnisky, March 2015.

viewingdirs=rotations(:,3,:); % Viewing direction is the third column of 
                              % each rotation matrix.
viewingdirs=squeeze(viewingdirs);
thetas=atan2d(viewingdirs(2,:),viewingdirs(1,:));
phis=acosd(viewingdirs(3,:));
%rhos=sqrt(1-viewingdirs(3,:).^2);

h1=figure;
%scatter(thetas(:),phis(:));
% The first parameters to polarscatter is indicating the angle in radians
% of each point, and the second is indicating the radius for each point. We
% use the polar angle as the radius.
% Since projection by an angle great than phi > pi/2 is the same as
% projection by projection by 180-phi plus in-plane rotation, we plot fold
% all phis into the range [0,pi/2].

folded_phis=phis;
idx=find(folded_phis>90);
folded_phis(idx)=180-folded_phis(idx);
assert(all(folded_phis<=90))
thetas(idx)=mod(thetas(idx)+180,180);

% Colors in the maps will be propotional to the density.
den=ksden(viewingdirs,0.2);
idx=floor((den-min(den))./(max(den)-min(den))*255)+1;
clrmap=jet(256);

if matlabversion >=9
    %polarscatter(thetas/180*pi,folded_phis,8,folded_phis,'filled');
    polarscatter(thetas/180*pi,folded_phis,8,clrmap(idx),'filled');
    rlim([0 90])
    rticks(0:10:90)
else
    %polar(thetas/180*pi,folded_phis,'o');
    polar(thetas/180*pi,clrmap(idx),'o');
end
clrbar=colorbar;
clrbar.Ticks=[0,0.5,1];
clrbar.TickLabels={'low','med','high'};

h2=figure;
scatter3(viewingdirs(1,:),viewingdirs(2,:),viewingdirs(3,:),5,clrmap(idx),'filled');
axis equal