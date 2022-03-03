% This script demonstrates how to match the Euler angles estimated by
% relation (and saved in the STAR files) and rotation matrices used by
% ASPIRE. It also demosntrates the relation between relion_project and
% cryo_project, and the relation between the estimated rotation and the
% Euler angles of relion.
% Finally, it shows how to fine for each estimated projetion its most
% similar ground truth image.
%
% Yoel Shkolnisky, March 2022.
 
star_fname='/data/yoelsh/datasets/10028/temp/Refine3D/job002/run_data.star';

% We use a STAR file from a refinement process since its contains the
% orientation paramteres for each of the images.
stardata_refined=readSTAR(star_fname);

% Read rotation matrices from the STAR file
N=5;
rots=zeros(3,3,N);
rots_t=zeros(3,3,N);
for k=1:N
    rec=stardata_refined(2).data{k};
    rots(:,:,k)=eul2rotm([rec.rlnAngleRot/180*pi,rec.rlnAngleTilt/180*pi,rec.rlnAnglePsi/180*pi],'ZYZ');
    rots_t(:,:,k)=(rots(:,:,k)).';
end

%% Read the volume, and project it using our projection code
vol_fname='/data/yoelsh/datasets/10028/temp/PostProcess/job004/postprocess.mrc';
vol=ReadMRC(vol_fname);
projs_aspire=cryo_project(vol,rots_t,129);  % Note that we use rots_t 
        % (transpose of what is written in the STAR file, since our code
        % internaly transposes the given rotations.
viewstack(projs_aspire,5,5)

%% Project the same volume using relion 
delete('proj.spi');
for k=1:N
    rec=stardata_refined(2).data{k};
    relion_command=sprintf('/opt/relion/bin/relion_project --i %s --rot %d --tilt %d --psi %d',...
        vol_fname,rec.rlnAngleRot,rec.rlnAngleTilt,rec.rlnAnglePsi);
    system(relion_command);
end
system('/opt/EMAN2/bin/e2proc2d.py proj.spi proj.mrcs');
projs_relion=ReadMRC('proj.mrcs');
delete('proj.mrcs');
viewstack(projs_relion,5,5) % Verify that projs_relion looks the same as 
    % projs_aspire. This implies that we know how to generate the same
    % projections in aspire and relion.

%% Now generate more projections, estimate their rotations, and compare to 
%% the ground truth  

% This is the same code as above.
N=500;
rots=zeros(3,3,N);
rots_t=zeros(3,3,N);
for k=1:N
    rec=stardata_refined(2).data{k};
    rots(:,:,k)=eul2rotm([rec.rlnAngleRot/180*pi,rec.rlnAngleTilt/180*pi,rec.rlnAnglePsi/180*pi],'ZYZ');
    rots_t(:,:,k)=(rots(:,:,k)).';
end
vol_fname='/data/yoelsh/datasets/10028/temp/PostProcess/job004/postprocess.mrc';
vol=ReadMRC(vol_fname);
projs_aspire=cryo_project(vol,rots_t,129);  % Note that we use rots_t 
        % (transpose of what is written in the STAR file, since our code
        % internaly transposes the given rotations.

% For consistency between the input projections and our reconstruction
% code, we need to first permute the dimensions of the images
projs_aspire=permute(projs_aspire,[2 1 3]);
        
projs_fname=tempmrcsname;
vol_fname=tempmrcname;
WriteMRC(projs_aspire,1.34*360/size(projs_aspire,1),projs_fname);
cryo_abinitio_C1_sync3n(projs_fname,vol_fname,'params.mat',1,1);
cryo_compare_volumes(vol_fname,2660,0.5,1.34*360/size(projs_aspire,1),1); %Sanity check
params=load('params.mat');
[regrot,mse,diff,O,flag]=register_rotations(params.rotations,rots);
fprintf('MSE of rotations estimation error is %7.3e\n',mse); % That is, we 
    % recovered the variable rots above generated from relion data. To plot
    % the view directions, we need generate rots from all data and plot it.
delete(projs_fname);
delete(vol_fname);


%% Plot view driretions for the entire dataset
% This is the same code as above.
N=length(stardata_refined(2).data);
rots=zeros(3,3,N);
rots_t=zeros(3,3,N);
for k=1:N
    rec=stardata_refined(2).data{k};
    rots(:,:,k)=eul2rotm([rec.rlnAngleRot/180*pi,rec.rlnAngleTilt/180*pi,rec.rlnAnglePsi/180*pi],'ZYZ');
    rots_t(:,:,k)=(rots(:,:,k)).';
end
cryo_plot_viewing_directions(rots);


%% Comapre denoised images to ground truth


% Generate clean reference projection from which we will choose the best
% matching projection for our image.
N=10000;
rots=rand_rots(N);
mapfile=cryo_fetch_emdID(8511);
volref=ReadMRC(mapfile);
volref=cryo_downsample(volref,129);
projs_ref=cryo_project(volref,rots);
projs_ref_normalized = zeros(size(projs_ref));
for k=1:N
    im=projs_ref(:,:,k);
    im=im-mean(im(:));
    projs_ref_normalized(:,:,k)=im./norm(im(:));
end

% Load denoised images
data=load('/data/Snir/ManifoldDenoising/Results/8511_simulation/20220220_8511_simulation_9000_SNR1.28.mat');

image_idx=[1,2,3];

g_projs_ref_normalized=gpuArray(single(projs_ref_normalized));
cmax=gpuArray(single(zeros(N,1)));

for k=1:numel(image_idx)
    idx=image_idx(k);
    % Noisy
    subplot(1,3,1);
    imagesc(data.clean_images(:,:,idx));
    colormap(gray);
    axis image;
    axis off;
    
    % Manifold denoising
    subplot(1,3,2);
    imagesc(data.denoised_projs_man_ctf(:,:,idx));
    colormap(gray);
    axis image;
    axis off;
    
    % Find the most similar image from the reference images
    im=data.clean_images(:,:,idx);
    %im=data.denoised_projs_man_ctf(:,:,idx);
    im=im-mean(im(:));
    im=im./norm(im(:));
    
    g_im=gpuArray(single(im));
    tic;
    
    parfor jj=1:N
        c=xcorr2(g_im,g_projs_ref_normalized(:,:,jj));
        cmax(jj)=max(c(:));
    end
    toc
    [bestcorr,bestidx]=max(cmax);
    subplot(1,3,3);
    imagesc(GaussFilt(projs_ref(:,:,bestidx),1));
    colormap(gray);
    axis image;
    axis off;
    disp(bestcorr)
    pause
end
