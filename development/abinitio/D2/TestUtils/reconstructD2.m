function [vol,dxD2] = reconstructD2(projs,rots,n_r,n_theta,max_shift,shift_step,dxD2_in,...
    shifts_2d_ref)

tic
disp('Reconstructing volume...');
[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
pfD2 = cat(3,pf,pf,pf,pf); % Replicate the images
%pfD2=pf;

if nargin>7
    shifts_ref_D2=repmat(shifts_2d_ref,4,1);
    %shifts_ref_D2=shifts_2d_ref;
else
    shifts_ref_D2=[];
end

gx = [1 0 0 ; 0 -1 0; 0 0 -1];
gy = [-1 0 0 ; 0 1 0; 0 0 -1];
gz = [-1 0 0 ; 0 -1 0; 0 0 1];

nImages = size(rots,3);

RsD2 = zeros(3,3,4*nImages);

for k=1:nImages
    rot = rots(:,:,k);
    RsD2(:,:,k)           = rot;
    RsD2(:,:,nImages+k)   = gx*rot;
    RsD2(:,:,2*nImages+k) = gy*rot;
    RsD2(:,:,3*nImages+k) = gz*rot;
end

if ~exist('dxD2_in','var') || (exist('dxD2_in','var') && isempty(dxD2_in))
    %[dxD2,~] = cryo_estimate_shifts(pf,rots,ceil(2*sqrt(2)*max_shift),shift_step,1,shifts_ref_D2,1);
    [dxD2,~] = ...
        cryo_estimate_shifts(pfD2,RsD2,ceil(2*sqrt(2)*max_shift),shift_step,1,shifts_ref_D2,1);
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
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsD2,RsD2,-dxD2, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);
toc