function err = measure_deviation_Cn(vol,n_symm)
   
errs = zeros(1,n_symm);
for i=1:n_symm-1
    vol_rot = fastrotate3z(vol,360*i/n_symm);
    errs(i) = norm(vol(:)-vol_rot(:))/norm(vol(:));
end
err = max(errs);

end