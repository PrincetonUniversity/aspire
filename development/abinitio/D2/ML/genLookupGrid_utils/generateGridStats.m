
function [grid_stats,pol_rep]=generateGridStats(rots_grid)
%% Generate rots_grid statistics for debuging
% latitude line angles and longitude values with respect to all symmetry 
% axes, angles between projection directions and 2 angles of spherical
% coordinates representation. 

% Project points to planes xy,xz,yz and calculate polar representation
% azimuth angle is the angle of latitude line
sphere_grid=squeeze(rots_grid(:,3,:));
plane_proj=projectToAllPlanes(sphere_grid);
azimuth_ang=[squeeze(atan2(sphere_grid(2,:),sphere_grid(1,:)))',...
            squeeze(atan2(sphere_grid(1,:),sphere_grid(3,:)))',...
            squeeze(atan2(sphere_grid(3,:),sphere_grid(2,:)))'];
pol_ang=squeeze(0.5*pi-acos(dot(plane_proj,repmat(sphere_grid,1,1,3),1)));
proj_dir_angles_z=acos(squeeze(rots_grid(:,3,:))'*squeeze(rots_grid(:,3,:)));
proj_dir_angles_y=acos(squeeze(rots_grid(:,2,:))'*squeeze(rots_grid(:,2,:)));
proj_dir_angles_x=acos(squeeze(rots_grid(:,1,:))'*squeeze(rots_grid(:,1,:)));

pol_rep=[azimuth_ang(:,1),pol_ang(:,1)];

azimuth_ang=azimuth_ang*180/pi;
pol_ang=pol_ang*180/pi;
proj_dir_angles_z=proj_dir_angles_z*180/pi;
proj_dir_angles_y=proj_dir_angles_y*180/pi;
proj_dir_angles_x=proj_dir_angles_x*180/pi;


%Containers for fields in statistics table
nrot=size(rots_grid,3);
npairs=nchoosek(nrot,2);
azimuth_x_i=zeros(npairs,1);
azimuth_y_i=zeros(npairs,1);
azimuth_z_i=zeros(npairs,1);
pol_x_i=zeros(npairs,1);
pol_y_i=zeros(npairs,1);
pol_z_i=zeros(npairs,1);
azimuth_x_j=zeros(npairs,1);
azimuth_y_j=zeros(npairs,1);
azimuth_z_j=zeros(npairs,1);
pol_x_j=zeros(npairs,1);
pol_y_j=zeros(npairs,1);
pol_z_j=zeros(npairs,1);

proj_angle_ij1=zeros(npairs,1);
proj_angle_ij2=zeros(npairs,1);
proj_angle_ij3=zeros(npairs,1);

azimuth_x_diff=zeros(npairs,1);
azimuth_y_diff=zeros(npairs,1);
azimuth_z_diff=zeros(npairs,1);
pol_x_diff=zeros(npairs,1);
pol_y_diff=zeros(npairs,1);
pol_z_diff=zeros(npairs,1);

%Fill in statistics
idx=0;
h=waitbar(0,'Generating grid statistics...');
for i=1:nrot-1
    waitbar(i/nrot-1);
    azimuth_x_i(idx+1:idx+nrot-i)=azimuth_ang(i,3);
    azimuth_y_i(idx+1:idx+nrot-i)=azimuth_ang(i,2);
    azimuth_z_i(idx+1:idx+nrot-i)=azimuth_ang(i,1);
    pol_x_i(idx+1:idx+nrot-i)=pol_ang(i,3);
    pol_y_i(idx+1:idx+nrot-i)=pol_ang(i,2);
    pol_z_i(idx+1:idx+nrot-i)=pol_ang(i,1);
    azimuth_x_j(idx+1:idx+nrot-i)=azimuth_ang(i+1:end,3);
    azimuth_y_j(idx+1:idx+nrot-i)=azimuth_ang(i+1:end,2);
    azimuth_z_j(idx+1:idx+nrot-i)=azimuth_ang(i+1:end,1);
    pol_x_j(idx+1:idx+nrot-i)=pol_ang(i+1:end,3);
    pol_y_j(idx+1:idx+nrot-i)=pol_ang(i+1:end,2);
    pol_z_j(idx+1:idx+nrot-i)=pol_ang(i+1:end,1);
    
    proj_angle_ij1(idx+1:idx+nrot-i)=proj_dir_angles_x(i,1+i:end);
    proj_angle_ij2(idx+1:idx+nrot-i)=proj_dir_angles_y(i,1+i:end);
    proj_angle_ij3(idx+1:idx+nrot-i)=proj_dir_angles_z(i,1+i:end);
    
    azimuth_x_diff(idx+1:idx+nrot-i)=abs(azimuth_ang(i,3)-azimuth_ang(i+1:end,3));
    azimuth_y_diff(idx+1:idx+nrot-i)=abs(azimuth_ang(i,2)-azimuth_ang(i+1:end,2));
    azimuth_z_diff(idx+1:idx+nrot-i)=abs(azimuth_ang(i,1)-azimuth_ang(i+1:end,1));
    pol_x_diff(idx+1:idx+nrot-i)=abs(pol_ang(i,3)-pol_ang(i+1:end,3));
    pol_y_diff(idx+1:idx+nrot-i)=abs(pol_ang(i,2)-pol_ang(i+1:end,2));
    pol_z_diff(idx+1:idx+nrot-i)=abs(pol_ang(i,1)-pol_ang(i+1:end,1));
    idx=idx+nrot-i;
end

ij_sub_idx=lin2sub_map(nrot);
row_names=char(strcat(repmat('(',npairs,1),string(ij_sub_idx(:,1)),...
    repmat(',',npairs,1),string(ij_sub_idx(:,2)),repmat(')',npairs,1)));
grid_stats=table(proj_angle_ij3,proj_angle_ij2,proj_angle_ij1,...
    azimuth_z_diff,azimuth_y_diff,azimuth_x_diff,pol_z_diff,pol_y_diff,...
    pol_x_diff,azimuth_z_i,azimuth_z_j,pol_z_i,pol_z_j,...
    azimuth_y_i,azimuth_y_j,pol_y_i,pol_y_j,azimuth_x_i,azimuth_x_j,...
    pol_x_i,pol_x_j,'RowNames',num2cell(row_names,2));
close(h);
end