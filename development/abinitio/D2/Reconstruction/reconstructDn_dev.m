function [vol,dxD2,corrs,Icl] = reconstructDn_dev(projs,rots,n_r,n_theta,max_shift,shift_step,dxD2_in,...
    shifts_2d_ref)

tic
disp('Reconstructing volume...');
[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
%pfD2 = cat(3,pf,pf,pf,pf); % Replicate the images
%pfD2=pf;

if exist('shifts_2d_ref','var')
    shifts_ref_D2=shifts_2d_ref;
    %shifts_ref_D2=shifts_2d_ref;
else
    shifts_ref_D2=[];
end

if ~exist('dxD2_in','var') || (exist('dxD2_in','var') && isempty(dxD2_in))
    %[dxD2,~] = cryo_estimate_shifts(pf,rots,ceil(2*sqrt(2)*max_shift),shift_step,1,shifts_ref_D2,1);
    [dxD2,RsD2,corrs,Icl] = ...
        cryo_estimate_shifts_Dn_dev(pf,rots,ceil(2*sqrt(2)*max_shift),shift_step,1,shifts_ref_D2,0);
    %[dxD2,~] = my_cryo_estimate_shifts(pfD2,RsD2,ceil(2*sqrt(2)*max_shift),shift_step,10000,shifts_ref_D2,1);
    %dxD2=-dxD2;
else
    dxD2=dxD2_in;
end
if max_shift==0
    dxD2=[];
end

%DEBUG
%dxD2=repmat(shifts_2d_ref,4,1);
%dxD2=repmat(dxD2,4,1);

n = size(projs,1);
projsD2 = cat(3,projs,projs,projs,projs);
dxD2=repmat(dxD2,4,1);
RsD2=cat(4,RsD2(:,:,:,1),RsD2(:,:,:,2),RsD2(:,:,:,3),RsD2(:,:,:,4));
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsD2,RsD2,-dxD2, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);
toc