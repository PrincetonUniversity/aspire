function vol = reconstruct_c3_c4(n_symm,projs,rots,n_r,n_theta,max_shift,shift_step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!!!This function is obsolete!!!!!!!!. Should replace it with reconstruct_cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_symm ~= 3 && n_symm ~= 4
    error('n_symm may be either 3 or 4');
end

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end

if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_r','var')
    n_r = ceil(size(projs,1)*0.5);
end

[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
if n_symm == 3
    pfC34 = cat(3,pf,pf,pf); % Replicate the images
else
    pfC34 = cat(3,pf,pf,pf,pf); % Replicate the images
end


g = [cosd(360/n_symm) -sind(360/n_symm) 0; ...
     sind(360/n_symm)  cosd(360/n_symm) 0; ...
     0                 0  1]; % rotation matrix of 120 or 90 degress around z-axis


nImages = size(rots,3);

RsC34 = zeros(3,3,n_symm*nImages);

for k=1:nImages
    rot = rots(:,:,k);
    RsC34(:,:,k)           = rot;
    RsC34(:,:,nImages+k)   = g*rot;
    RsC34(:,:,2*nImages+k) = g*g*rot;
    if n_symm == 4
        RsC34(:,:,3*nImages+k) = g*g*g*rot;
    end
end

[dxC34,~] = cryo_estimate_shifts(pfC34,RsC34,max_shift,shift_step,10000,[],0);

n = size(projs,1);
if n_symm == 3
    projsC34 = cat(3,projs,projs,projs);
else
    projsC34 = cat(3,projs,projs,projs,projs);
end

[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projsC34,RsC34,-dxC34, 1e-6, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);