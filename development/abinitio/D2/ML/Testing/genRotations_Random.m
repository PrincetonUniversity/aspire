

%Input: N is size of grid, best choose N=k^2
function [debug_data]=genRotations_Random(N,eq_filter_angle,s)
%% Generate random grid on the sphere
%rng(s);
q=qrand(N);
rots_grid=zeros(3,3,N);
for i=1:N
    rots_grid(:,:,i)=q_to_rot(q(:,i));
end
[sphere_grid]=squeeze(rots_grid(:,3,:))';
% octant1_idx=(sphere_grid(:,1)>0).*(sphere_grid(:,2)>0).*(sphere_grid(:,3)>0);
% sphere_grid=sphere_grid(octant1_idx==1,:);
% rots_grid=rots_grid(:,:,octant1_idx==1);
nrot=size(sphere_grid,1);

%% Filter points on equators
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
filtered_idx=~sum(angular_dists>eq_min_dist,2)>0;
sphere_grid=sphere_grid(:,filtered_idx);
rots_grid=rots_grid(:,:,filtered_idx);
q=q(:,filtered_idx);
nrot=size(sphere_grid,2);
clearvars filtered_idx eq_min_dist angular_dists %proj_xy proj_xz proj yz

%% Generate rots_grid statistics for debuging
% latitude line angles and longitude values with respect to all symmetry 
% axes, angles between projection directions and 2 angles of spherical
% coordinates representation. 
[grid_stats,pol_rep]=generateGridStats(rots_grid);

%% Finalize
debug_data=struct('rots_grid',rots_grid,'grid_stats',grid_stats,'nrot',nrot,...
    'eq_filter_angle',eq_filter_angle,'pol_rep',pol_rep,'q',q);
end