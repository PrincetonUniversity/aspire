function [projections]=shift_images(projections, shifts)
%This function use FFT to shift images
%Input:     projections: LxLxK, K is the number of images
%           shifts: Kx2 array, K is the number of images
%Output:    projections: LxLxP
%Zhizhen Zhao 09/01/2012

%shift the images in cartesian grid number of projections is the same as
%the length of shifts

K=size(projections, 3);%number of images
n=size(projections, 1);%size of the image
range=-fix(n/2):fix(n/2);
[omega_x,omega_y]=ndgrid(range,range);
omega_x=-2*pi.*omega_x/n; omega_y=-2*pi.*omega_y/n;

%shift each image
for k=1:K
    p=projections(:,:,k);
    pf=fftshift(fft2(ifftshift(p)));
    phase_x=omega_x.*shifts(k,1);
    phase_y=omega_y.*shifts(k,2);
    pf=pf.*exp(sqrt(-1)*(phase_x+phase_y));
    p2=fftshift(ifft2(ifftshift(pf)));
    projections(:,:,k)=real(p2);
end

end