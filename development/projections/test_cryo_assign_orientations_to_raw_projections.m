% XXX Under construction!
% XXX This file should eventually test the function cryo_assign_orientations_to_raw_projections.
% XXX The current version is some draft used for development.
%
% Yoel Shkolnisky, July 2016.

clear

% Load data
alignment_data=load('/home/yoel/tmp/cleansim/averages_info_nn10_group1.mat');
averaging_data=load('/home/yoel/tmp/cleansim/abinitio_info_nn10_nm1000_group1.mat');

averaging_data.nnavg=10;
% Compute reference rotations
% Load reference rotations
load /home/yoel/tmp/cleansim/cleandata.mat q shifts
qref=q;
Rref=zeros(3,3,size(qref,2));
for j=1:size(Rref,3)
    R=q_to_rot(qref(:,j));
    Rref(:,:,j)=R.';
end

[estR,est_shifts,rawimageindices]=cryo_assign_orientations_to_raw_projections(...
    alignment_data,averaging_data,Rref,shifts);

% [estR,est_shifts,rawimageindices]=cryo_assign_orientations_to_raw_projections(...
%     alignment_data,averaging_data);

projections=ReadMRC('/home/yoel/tmp/cleansim/phaseflipped_cropped_downsampled_prewhitened_group1.mrc');
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projections(:,:,rawimageindices),estR,-est_shifts, 1e-6, 100, zeros(129,129,129));
[ v2, ~, ~ ,~, ~, ~] = recon3d_firm( projections(:,:,rawimageindices),Rref(:,:,rawimageindices),-shifts(rawimageindices,:), 1e-6, 100, zeros(129,129,129));
[estR,estdx,vol2aligned,reflect]=cryo_align_densities(v1,v2,1,1);

plotFSC(v1,vol2aligned,0.143,1)
