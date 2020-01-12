
%% Input:  sphere_grid== grid of points on 3D sphere
%          eq_filter_angle== angular distance from equator to be marked as
%          and equator point
%          n= n of Dn. 
%  Output: eq_idx= indices of points on sphere which distance from one of
%  the equators is < eq_filter angle

function [eq_idx,eq_class,angular_dists]=markEquatorsDn2(sphere_grid,eq_filter_angles,n) 
%% project each vector onto xy,xz,yz planes and measure angular distance
nrot=size(sphere_grid,2);
angular_dists=zeros(nrot,3);
proj_xy=sphere_grid;
proj_xy(3,:)=0;
norms_xy=sqrt(sum(proj_xy(1:2,:).^2,1));
proj_xy=proj_xy./repmat(norms_xy,3,1);
angular_dists(:,1)=dot(sphere_grid,proj_xy,1);

b2=[0,0,1];
coeff2=b2*sphere_grid;

if mod(n,2)>0
theta=(n-1)/2*pi/n;
b1=[-sin(theta),cos(theta),0];
coeff1=b1*sphere_grid;
proj_xz=repmat(coeff1,3,1).*b1'+repmat(coeff2,3,1).*b2';
norms_xz=sqrt(sum(proj_xz(1:3,:).^2,1));
proj_xz=proj_xz./repmat(norms_xz,3,1);
angular_dists(:,2)=dot(sphere_grid,proj_xz,1);

theta=(n+1)/2*pi/n;
b1=[-sin(theta),cos(theta),0];
coeff1=b1*sphere_grid;
proj_yz=repmat(coeff1,3,1).*b1'+repmat(coeff2,3,1).*b2';
norms_yz=sqrt(sum(proj_yz(1:3,:).^2,1));
proj_yz=proj_yz./repmat(norms_yz,3,1);
angular_dists(:,3)=dot(sphere_grid,proj_yz,1);

else

proj_xz=sphere_grid;
proj_xz(2,:)=0;
norms_xz=sqrt(sum(proj_xz([1,3],:).^2,1));
proj_xz=proj_xz./repmat(norms_xz,3,1);
angular_dists(:,2)=dot(sphere_grid,proj_xz,1);

theta=pi/n;
b1=[cos(theta),sin(theta),0];
coeff1=b1*sphere_grid;
proj_yz=repmat(coeff1,3,1).*b1'+repmat(coeff2,3,1).*b2';
norms_yz=sqrt(sum(proj_yz(1:3,:).^2,1));
proj_yz=proj_yz./repmat(norms_yz,3,1);
angular_dists(:,3)=dot(sphere_grid,proj_yz,1);

end

%Mark points close to equators
eq_min_dist=cos(eq_filter_angles(1)*pi/180);
eq_idx_z=sum(angular_dists(:,1)>eq_min_dist,2)>0;
eq_min_dist=cos(eq_filter_angles(2)*pi/180);
eq_idx_y=sum(angular_dists(:,2)>eq_min_dist,2)>0;
eq_min_dist=cos(eq_filter_angles(3)*pi/180);
eq_idx_x=sum(angular_dists(:,3)>eq_min_dist,2)>0;
eq_idx=[eq_idx_z,eq_idx_y,eq_idx_x];

%Classify equators
% 1 -> z equator
% 2 -> y equator
% 3 -> x equator
% 4 -> z top view, i.e. both x and y equator
% 5 -> y top view, i.e. both x and z equator
% 6 -> x top view, i.e. both y and z equator
eq_class=zeros(nrot,1);
%nEq=sum(angular_dists>eq_min_dist,2);
nEq=sum(eq_idx,2);
topView_idx=nEq>1;
[~,topView_class]=min(eq_idx(topView_idx,:),[],2);
eq_class(topView_idx)=topView_class+3;
non_topView_idx=(nEq==1);
%[~,non_topView_class]=max(angular_dists(non_topView_idx,:)>eq_min_dist,[],2);
[~,non_topView_class]=max(eq_idx(non_topView_idx,:),[],2);
eq_class(non_topView_idx)=non_topView_class;
eq_idx=max(eq_idx,[],2);





