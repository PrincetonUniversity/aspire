function vol_out = make_vol_c12(vol_in)

vol_out = vol_in;
for i=1:11
    vol_out = vol_out + fastrotate3z(vol_in,i*360/12);
end
vol_out = vol_out/12;

end