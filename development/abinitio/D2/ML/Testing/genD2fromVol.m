
function [D2_vol]=genD2fromVol(vol)

vol_z=flip(flip(vol,1),2);
vol_y=flip(flip(vol,1),3);
vol_x=flip(flip(vol,2),3);
D2_vol=vol+vol_x+vol_y+vol_z;