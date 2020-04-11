function [D2_vol]=genD2fromVol(vol)
% GEND2VOLFROMVOL   Enforce volume to be D2
%
% [D2_vol]=genD2fromVol(vol)
%   Enforce vol to be D2 symmetric by averaging it with its copies rotated
%   by all D2 group elements.
%
% Eitan Rosen, January 2020.

vol_z=flip(flip(vol,1),2);
vol_y=flip(flip(vol,1),3);
vol_x=flip(flip(vol,2),3);
D2_vol=(vol+vol_x+vol_y+vol_z)/4;