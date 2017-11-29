function err = remove_assertVolumeIsC12(vol)
   
errs = zeros(1,12);
for i=1:11
    vol_rot = fastrotate3z(vol,360*i/12);
    errs(i) = norm(vol(:)-vol_rot(:))/norm(vol(:));
end
err = max(errs);
log_message('deviation of volume from c12 symmetry is %.4f',err);
end