projs_mrc='./datasets/10081/averages_nn50_group1.mrc';
vol_mrc='vol_10081_nn50_grp1_5000.mrc';
map_out_step1='c4_10081_map_out_1_group1.mrc';
map_out_step2='c4_10081_map_out_2_group1.mrc';
mat_out='c4_10081_map_out_group1.mat';
cryo_refine_map(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out,25);

projs_mrc='./datasets/10081/averages_nn50_group2.mrc';
vol_mrc='vol_10081_nn50_grp2_5000.mrc';
map_out_step1='c4_10081_map_out_1_group2.mrc';
map_out_step2='c4_10081_map_out_2_group2.mrc';
mat_out='c4_10081_map_out_group2.mat';
cryo_refine_map(projs_mrc,vol_mrc,map_out_step1,map_out_step2,25,mat_out);

% vol1a=ReadMRC('c4_10081_map_out_1_group1.mrc');
% vol2a=ReadMRC('c4_10081_map_out_1_group2.mrc');
% [Rest,estdx,vol2aaligned]=cryo_align_densities(vol1a,vol2a,2.5698,1,[],0,50);
% plotFSC(vol1a,vol2aaligned,0.143,2.5698);
% WriteMRC(vol2aaligned,1,'c4_10081_map_out_1_group2_aligned_to_group1.mrc');

vol1b=ReadMRC('c4_10081_map_out_2_group1.mrc');
vol2b=ReadMRC('c4_10081_map_out_2_group2.mrc');
[Rest,estdx,vol2baligned]=cryo_align_densities(vol1b,vol2b,2.5698,1,[],0,50);
plotFSC(vol1b,vol2baligned,0.143,2.5698);
WriteMRC(vol2baligned,1,'c4_10081_map_out_2_group2_aligned_to_group1.mrc');