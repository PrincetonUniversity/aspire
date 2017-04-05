% vol_g1='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group1.mrc';
% projs_g1='/home/yoel/Desktop/tmp/80S_89/averages_nn10_group1.mrc';
% p1=ReadMRC(projs_g1);
% projs_g1=tempmrcname;
% WriteMRC(p1(:,:,1:5000),1,projs_g1);
% clear p1
% 
% vol_g2='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group2.mrc';
% projs_g2='/home/yoel/Desktop/tmp/80S_89/averages_nn10_group2.mrc';
% p2=ReadMRC(projs_g2);
% projs_g2=tempmrcname;
% WriteMRC(p2(:,:,1:5000),1,projs_g2);
% clear p2

% vol_g1='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group1.mrc';
% projs_g1='/home/yoel/Desktop/tmp/80S_129/averages_nn10_group1.mrc';
% p1=ReadMRC(projs_g1);
% projs_g1=tempmrcname;
% WriteMRC(p1(:,:,1:10000),1,projs_g1);
% clear p1
% 
% vol_g2='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group2.mrc';
% projs_g2='/home/yoel/Desktop/tmp/80S_129/averages_nn10_group2.mrc';
% p2=ReadMRC(projs_g2);
% projs_g2=tempmrcname;
% WriteMRC(p2(:,:,1:10000),1,projs_g2);
% clear p2

vol_g1='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group1.mrc';
projs_g1='/home/yoel/Desktop/tmp/80S_129/phaseflipped_downsampled_prewhitened_group1.mrc';
p1=ReadMRC(projs_g1);
projs_g1=tempmrcname;
WriteMRC(p1(:,:,1:end),1,projs_g1);
clear p1

vol_g2='/home/yoel/Desktop/tmp/80S_89/vol_nn10_nm2000_group2.mrc';
projs_g2='/home/yoel/Desktop/tmp/80S_129/phaseflipped_downsampled_prewhitened_group2.mrc';
p2=ReadMRC(projs_g2);
projs_g2=tempmrcname;
WriteMRC(p2(:,:,1:end),1,projs_g2);
clear p2

% Align densities and plot FSC
pixA=1.43*359/89;
v1=ReadMRC(vol_g1); v2=ReadMRC(vol_g2);
[~,~,v2alignedtov1]=cryo_align_densities(v1,v2,pixA,1);
plotFSC(v1,v2alignedtov1,0.143,pixA);

% Refine g1
map_out_step1_g1=fullfile(tempmrcdir,'map_out_step1_g1');
map_out_step2_g1=fullfile(tempmrcdir,'map_out_step2_g1');
mat_out_g1=fullfile(tempmrcdir,'mat_out_g1');
cryo_refine_map_iter(projs_g1,vol_g1,map_out_step1_g1,map_out_step2_g1,mat_out_g1)

% Refine g2
map_out_step1_g2=fullfile(tempmrcdir,'map_out_step1_g2');
map_out_step2_g2=fullfile(tempmrcdir,'map_out_step2_g2');
mat_out_g2=fullfile(tempmrcdir,'mat_out_g2');
cryo_refine_map_iter(projs_g2,vol_g2,map_out_step1_g2,map_out_step2_g2,mat_out_g2)