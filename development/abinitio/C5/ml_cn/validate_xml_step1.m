
function validate_xml_step1(xmlname_full)
% Read workflow file
tree = xmltree(xmlname_full);
workflow = convert(tree);
% Validate struct

assertfield(workflow,'info','name');
assertfield(workflow,'info','description');
assertfield(workflow,'info','created');
assertfield(workflow,'info','working_dir');
assertfield(workflow,'info','logfile');

assertfield(workflow,'algo','n_symm');
assertfield(workflow,'algo','n_images');
assertfield(workflow,'algo','image_size');
assertfield(workflow,'algo','snr');
assertfield(workflow,'algo','n_r_perc');
assertfield(workflow,'algo','max_shift_perc');
assertfield(workflow,'algo','shift_step');
assertfield(workflow,'algo','do_handle_equators');
assertfield(workflow,'algo','inplane_rot_res');
assertfield(workflow,'algo','do_recons_vol');
assertfield(workflow,'algo','recon_mrc_fname');
assertfield(workflow,'algo','recon_mat_fname');

assertfield(workflow,'cache','name');

end