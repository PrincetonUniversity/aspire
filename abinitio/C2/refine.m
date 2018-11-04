projs_mrc = '/home/yoel/scratch/10038_refspick/eman/particles__ctf_flip.mrcs';
vol_mrc   = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C2/10038_refspick_nn50_ims500.mrc';
map_out_step1 =  'map_out_step1_10038_refspick.mrc';
map_out_step2 =  'map_out_step2_10038_refspick.mrc';
mat_out = 'map_out_10038_refspick_maxiter10';
maxiter  = 10;

cryo_refine_map_c2(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out,maxiter);