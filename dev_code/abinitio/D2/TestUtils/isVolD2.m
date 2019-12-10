
function [err]=isVolD2(vol)
%% Rotate in all directions
vol_z=flip(flip(vol,1),2);
vol_y=flip(flip(vol,1),3);
vol_x=flip(flip(vol,2),3);

%err=abs(vol-vol_z)+abs(vol-vol_y)+abs(vol-vol_x);
nvol=sqrt(sum(vol(:).^2));
nz=(vol-vol_z).^2;
ny=(vol-vol_y).^2;
nx=(vol-vol_x).^2;
err=[sqrt(sum(nz(:))),sqrt(sum(ny(:))),sqrt(sum(nx(:)))]/nvol;

