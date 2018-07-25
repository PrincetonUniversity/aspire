function save_vols(vol,vol_filename,n_symm)

vol_filt    = make_vol_cn(vol,n_symm,true);
vol_no_filt = make_vol_cn(vol,n_symm,false);

[folder, baseFileName, extension] = fileparts(vol_filename);
vol_filt_file_name = sprintf('%s_reinf_filt%s',fullfile(folder, baseFileName),extension);
vol_no_filt_file_name = sprintf('%s_reinf%s',fullfile(folder, baseFileName),extension);

err = measure_deviation_Cn(vol,n_symm);
log_message('deviation of abinitio volume from C%d symmetry is %.4f',n_symm,err);
err = measure_deviation_Cn(vol_filt,n_symm);
log_message('deviation of abinitio C%d-enforced filtered volume from C%d symmetry is %.4f',n_symm,n_symm,err);
err = measure_deviation_Cn(vol_no_filt,n_symm);
log_message('deviation of abinitio C%d-enforced non-filtered volume from C%d symmetry is %.4f',n_symm,n_symm,err);

log_message('Saving abinitio volume under: %s', vol_filename);
WriteMRC(vol,1,vol_filename);

log_message('Saving abinitio C%d-enforced filtered volume under: %s', n_symm,vol_filt_file_name);
WriteMRC(vol_filt,1,vol_filt_file_name);

log_message('Saving abinitio C%d-enforced non-filtered volume under: %s', n_symm,vol_no_filt_file_name);
WriteMRC(vol_no_filt,1,vol_no_filt_file_name);

end


% function save_vols(vol,vol_filename,n_symm)
% 
% assertVolumeIsCn(vol,n_symm);
% log_message('Saving abinitio volume under: %s', vol_filename);
% WriteMRC(vol,1,vol_filename);
% 
% vol_filt    = make_vol_cn(vol,n_symm,true);
% assertVolumeIsCn(vol_filt,n_symm);
% 
% vol_no_filt = make_vol_cn(vol,n_symm,false);
% assertVolumeIsCn(vol_no_filt,n_symm);
% 
% [folder, baseFileName, extension] = fileparts(vol_filename);
% 
% filt_file_name = sprintf('%s_reinf_filt%s',fullfile(folder, baseFileName),extension);
% WriteMRC(vol_filt,1,filt_file_name);
% 
% no_filt_file_name = sprintf('%s_reinf%s',fullfile(folder, baseFileName),extension);
% WriteMRC(vol_no_filt,1,no_filt_file_name);
% 
% end