% xmlname_full = cryo_workflow_abinitio_Cn_ml_start;
% 
% %% option A. Call using an xml
% cryo_workflow_abinitio_Cn_ml_execute(xmlname_full);

%% option B. Call using explicit variables
% Read workflow file
tree = xmltree(xmlname_full);
workflow = convert(tree);

% A call to impl with only mandatory variables
n_symm = str2double(workflow.algo.n_symm);
mrc_stack_file = workflow.info.mrc_stack_file;
recon_mrc_fname = workflow.algo.recon_mrc_fname;

cryo_workflow_abinitio_Cn_ml_execute2(n_symm,mrc_stack_file,recon_mrc_fname);

% A call to impl with all possible variables
cache_file_name = workflow.cache.name;
recon_mat_fname = workflow.algo.recon_mat_fname;
do_downsample = str2double(workflow.algo.do_downsample);
downsample_size = str2double(workflow.algo.downsample_size);
n_r_perc = str2double(workflow.algo.n_r_perc);
max_shift_perc = str2double(workflow.algo.max_shift_perc);
shift_step = str2double(workflow.algo.shift_step);
mask_radius_perc = str2double(workflow.algo.mask_radius_perc);
do_handle_equators = str2double(workflow.algo.do_handle_equators);
inplane_rot_res = str2double(workflow.algo.inplane_rot_res);

cryo_workflow_abinitio_Cn_ml_execute2(n_symm,mrc_stack_file,recon_mrc_fname,cache_file_name,recon_mat_fname,...
    do_downsample,downsample_size,n_r_perc,max_shift_perc,shift_step,mask_radius_perc,do_handle_equators,inplane_rot_res);