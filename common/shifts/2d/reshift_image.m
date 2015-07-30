function sim=reshift_image(im,s)
%
% Shift an image im by the vector s.
% s is a vector of two coordinates, whose first coordinate is the number of pixels to translate in the x direction, and whose second coordinate in the number of pixels to translate in the y direction.
%
% s DOES NOT need be a vector of integer, so translation of fractional pixels is allowed. im can be of odd or even side.
%
% Note: I don't know if s=[1 0 ] will translate left or right. Just check it. Same for s=[0 1].
%
% Example:   reshift_image2(im,[-0.5 0 ])
%
% This function replaces reshift_image.m which works only for images with odd side.
%
% Yoel Shkolnisky, August 2012.


if ~ismatrix(im)
    error('Input must be a 2D image');
end

if size(im,1)~=size(im,2)
    error('Image must be square');
end

n=size(im,1);
ll=fix(n/2);
freqrng=-ll:n-ll-1;
[omega_x,omega_y]=ndgrid(freqrng,freqrng);
omega_x=2*pi.*omega_x/n;
omega_y=2*pi.*omega_y/n;

phase_x=exp(sqrt(-1)*omega_x.*s(1));
phase_y=exp(sqrt(-1)*omega_y.*s(2));

% Force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod(n, 2) == 0
    phase_x(1,:)=real(phase_x(1,:));
    phase_y(:,1)=real(phase_y(:,1));
end 

phases=phase_x.*phase_y;
pim=fftshift(fft2(ifftshift(im)));
pim=pim.*phases;
sim=fftshift(ifft2(ifftshift(pim)));

if norm(imag(sim(:)))/norm(sim(:))>1.0e-8
    error('Large imaginary components');
end
sim=real(sim);

