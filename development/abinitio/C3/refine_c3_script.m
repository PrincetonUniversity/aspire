projs_mrc = '/home/yoel/Desktop/gabyf/output/averages_nn50_group2.mrc';
% vol_mrc   = 'c3Frank_nn50_grp2_ims5000.mrc';
vol_mrc = 'model_07_01_1p1.mrc';
map_out_step1 =  'map_out_step1_model_07_01_1p1.mrc';
map_out_step2 =  'map_out_step2_model_07_01_1p1.mrc';
mat_out = 'map_out_model_07_01_1p1.mrc';
maxiter  = 5; 
cryo_refine_map_c3(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out,maxiter);