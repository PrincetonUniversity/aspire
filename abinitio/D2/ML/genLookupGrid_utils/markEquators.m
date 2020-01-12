
%% Input:  sphere_grid== grid of points on 3D sphere
%          eq_filter_angle== angular distance from equator to be marked as
%          and equator point
%  Output: eq_idx= indices of points on sphere which distance from one of
%  the equators is < eq_filter angle

function [eq_idx,eq_class,angular_dists]=markEquators(sphere_grid,eq_filter_angle) 
%% project each vector onto xy,xz,yz planes and measure angular distance
nrot=size(sphere_grid,2);
angular_dists=zeros(nrot,3);
proj_xy=sphere_grid;
proj_xy(3,:)=0;
norms_xy=sqrt(sum(proj_xy(1:2,:).^2,1));
proj_xy=proj_xy./repmat(norms_xy,3,1);
angular_dists(:,1)=dot(sphere_grid,proj_xy,1);
proj_xz=sphere_grid;
proj_xz(2,:)=0;
norms_xz=sqrt(sum(proj_xz([1,3],:).^2,1));
proj_xz=proj_xz./repmat(norms_xz,3,1);
angular_dists(:,2)=dot(sphere_grid,proj_xz,1);
proj_yz=sphere_grid;
proj_yz(1,:)=0;
norms_yz=sqrt(sum(proj_yz(2:3,:).^2,1));
proj_yz=proj_yz./repmat(norms_yz,3,1);
angular_dists(:,3)=dot(sphere_grid,proj_yz,1);

%Mark points close to equators (within distance eq_filter_angle)
eq_min_dist=cos(eq_filter_angle*pi/180);
eq_idx=sum(angular_dists>eq_min_dist,2)>0;

%Classify equators
% 1 -> z equator
% 2 -> y equator
% 3 -> x equator
% 4 -> z top view, i.e. both x and y equator
% 5 -> y top view, i.e. both x and z equator
% 6 -> x top view, i.e. both y and z equator
eq_class=zeros(nrot,1);
nEq=sum(angular_dists>eq_min_dist,2);
topView_idx=nEq>1;
[~,topView_class]=min(angular_dists(topView_idx,:)>eq_min_dist,[],2);
eq_class(topView_idx)=topView_class+3;
non_topView_idx=(nEq==1);
[~,non_topView_class]=max(angular_dists(non_topView_idx,:)>eq_min_dist,[],2);
eq_class(non_topView_idx)=non_topView_class;





