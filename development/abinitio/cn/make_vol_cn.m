function vol_out = make_vol_cn(vol_in,n_symm,filt)

if filt
%     vol_in = cryo_downsample(vol_in,89);
%     vol_in = GaussFilt(vol_in,0.2);
    sz = size(vol_in);
    vol_in = cryo_mask_volume(vol_in,sz(1)*0.45,sz(1)*0.05);
end

vol_out = vol_in;
for i=1:n_symm-1
    vol_out = vol_out + fastrotate3z(vol_in,i*360/n_symm);
end
vol_out = vol_out/n_symm;

% if filt
%     vol_out = cryo_mask_volume(vol_out,25,5);
% end

end

    