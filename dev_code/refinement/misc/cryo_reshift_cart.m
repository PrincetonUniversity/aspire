function [ projections ] = cryo_reshift_cart( projections, shift )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%   projections are the non-centered images and shifts are their relative
%   shifts to the center of mass. This program reshift images to center
%   according to the value of the shifts.

N=ceil(size(projections, 1)/2); K=size(projections, 3);
[omega_x,omega_y]=ndgrid(-(N-1):N-1,-(N-1):N-1);
omega_x=-2*pi.*omega_x/(2*N-1); omega_y=-2*pi.*omega_y/(2*N-1);

for k=1:K
    p=projections(:,:,k);
    pf=cfftn(p);
    phase_x=omega_x.*shift(1);
    phase_y=omega_y.*shift(2);
    pf=pf.*exp(-sqrt(-1)*(phase_x+phase_y));
    p2=icfftn(pf);
    projections(:,:,k)=real(p2);
end

end

