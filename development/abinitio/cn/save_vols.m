function [vol_filename,vol_filt_file_name,vol_no_filt_file_name] = save_vols(vol,vol_filename,n_symm)

vol_filt    = make_vol_cn(vol,n_symm,true);
vol_no_filt = make_vol_cn(vol,n_symm,false);

[folder, baseFileName, extension] = fileparts(vol_filename);
vol_filt_file_name = sprintf('%s_reinf_filt%s',fullfile(folder, baseFileName),extension);
vol_no_filt_file_name = sprintf('%s_reinf%s',fullfile(folder, baseFileName),extension);

err = measure_deviation_cn(vol,n_symm);
log_message('deviation of abinitio volume from C%d symmetry is %.4f',n_symm,err);
err = measure_deviation_cn(vol_filt,n_symm);
log_message('deviation of abinitio C%d-enforced filtered volume from C%d symmetry is %.4f',n_symm,n_symm,err);
err = measure_deviation_cn(vol_no_filt,n_symm);
log_message('deviation of abinitio C%d-enforced non-filtered volume from C%d symmetry is %.4f',n_symm,n_symm,err);

log_message('Saving abinitio volume under: %s', vol_filename);
WriteMRC(vol,1,vol_filename);

log_message('Saving abinitio C%d-enforced filtered volume under: %s', n_symm,vol_filt_file_name);
WriteMRC(vol_filt,1,vol_filt_file_name);

log_message('Saving abinitio C%d-enforced non-filtered volume under: %s', n_symm,vol_no_filt_file_name);
WriteMRC(vol_no_filt,1,vol_no_filt_file_name);

end


function err = measure_deviation_cn(vol,n_symm)
   
errs = zeros(1,n_symm);
for i=1:n_symm-1
    vol_rot = fastrotate3z(vol,360*i/n_symm);
    errs(i) = norm(vol(:)-vol_rot(:))/norm(vol(:));
end
err = max(errs);

end