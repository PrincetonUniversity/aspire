function svol=reshift_vol(vol,s)
%
% Shift the volume given by im by the vector s using trigonometric
% interpolation. The volume im is of nxnxn, where n can be odi or even. The vector
% s\in \mathbb{R}^{3} contains the hshifts.
%
% Example: Shift the volume vol by 1 pixel in the x direction, 2 in the y
% direction, and 3 in the z direction
%
%       s = [1 2 3];
%       vols=shift_vol(vol,s);
%
% NOTE: I don't know if s=[0 0 1 ] shifts up or down, but this can be easily checked. Same issue for the other directions.
%
% Yoe Shkolnisky, August 2012.

if ndims(vol)~=3
    error('Input must be a 3D volume');
end

if (size(vol,1)~=size(vol,2)) || (size(vol,1)~=size(vol,3))
    error('All three dimension of the input must be equal');
end

n=size(vol,1);
ll=fix(n/2);
freqrng=-ll:n-ll-1;
[omega_x,omega_y,omega_z]=ndgrid(freqrng,freqrng,freqrng);
omega_x=2*pi.*omega_x/n; 
omega_y=2*pi.*omega_y/n;
omega_z=2*pi.*omega_z/n;

phase_x=exp(1i.*omega_x.*s(1));
phase_y=exp(1i.*omega_y.*s(2));
phase_z=exp(1i.*omega_z.*s(3));


% Force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod(n,2)==0
    phase_x(1,:,:)=real(phase_x(1,:,:));
    phase_y(:,1,:)=real(phase_y(:,1,:));
    phase_z(:,:,1)=real(phase_z(:,:,1));    
end

phases=phase_x.*phase_y.*phase_z;
pim=fftshift(fftn(ifftshift(vol)));
pim=pim.*phases;
svol=fftshift(ifftn(ifftshift(pim)));

if norm(imag(svol(:)))/norm(svol(:))>1.0e-8
    error('Large imaginary components');
end
svol=real(svol);
