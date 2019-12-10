
function [eq_idx,eq_class,angular_dists]=cryo_detect_equator_rots(rots,eq_filter_angle)

[sphere_grid]=squeeze(rots(:,3,:))';
nrot=size(sphere_grid,1);

%project each vector onto xy,xz,yz planes and measure angular distance
sphere_grid=sphere_grid'; 
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

%filter points close to equators
eq_min_dist=cos(eq_filter_angle*pi/180);
eq_idx=sum(angular_dists>eq_min_dist,2)>0;%find(sum(angular_dists>eq_min_dist,2)>0);
eq_class=angular_dists>eq_min_dist;
%[~,eq_class]=max(angular_dists(eq_idx,:),[],2);

