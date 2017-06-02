function cryo_refine_map_outofcore(projs_fname,vol_fname,map_out_step1,map_out_step2,mat_out,maxiter,tol)

if ~exist('maxiter','var')
    maxiter=5;
end

if ~exist('tol','var')
    tol=10;
end

% Strip extension from MRC filenames
[d,n,~]=fileparts(map_out_step2);
map_out_step2=fullfile(d,n);

cr_threshold=0.8;

log_message('Using maxiter=%d, tol=%4.2f degrees, correlation threshold=%d;',maxiter,tol,cr_threshold);
log_message('Iterations stop when maxiter reached or the rotations of the current iteration deviate from the rotations of the previous');
log_message('iteatrion by a less than tol degrees, or the correlation between volumes of two consecutive iterations is at least as specified.');

log_message('Loading initial volume %s',vol_fname);
vol=ReadMRC(vol_fname);
log_message('Volume loaded');
szvol=size(vol);

log_message('Applying global phase-flip to initial volume');
vol=cryo_globalphaseflip_vol(vol);
log_message('Volume phase-flipped');

log_message('Applying global phase-flip projections');
projs_fname_phaseflipped=tempmrcname;
cryo_globalphaseflip_outofcore(projs_fname,projs_fname_phaseflipped);
log_message('Projections phase-flipped');

projs=imagestackReader(projs_fname_phaseflipped,100);
szprojs=projs.dim;


log_message('Dimensions of volume %dx%dx%d]',szvol(1),szvol(2),szvol(3));
log_message('Dimensions of projections %dx%d (%d projections)',szprojs(1),szprojs(2),szprojs(3));
if szvol(1)~=szprojs(1)
    log_message('Volume and projections have differnt dimensions. Resampling initial volume to size %dx%dx%d',...
        szprojs(1),szprojs(1),szprojs(1));
    vol=cryo_downsample(vol,[szprojs(1),szprojs(1),szprojs(1)]);
    log_message('Initial volume resampled to %dx%dx%d',szprojs(1),szprojs(1),szprojs(1));
end
log_message('vol MD5 %s',MD5var(vol));

szvol=size(vol); % Dimensions of volume may have changed due to resampling
% Preprocess projections and referece volume.
log_message('Start preprocessing projections');
n=szvol(1);
projs_fname_normalized=tempmrcname;
cryo_mask_outofcore(projs_fname_phaseflipped,projs_fname_normalized,floor(0.45*n),floor(0.05*n)); % Mask projections
log_message('Preprocessing done');
delete(projs_fname_phaseflipped);

Nrefs=100;
log_message('Using Nrefs=%d reference projections for alignment',Nrefs);

% Set the angular resolution for common lines calculations. The resolution
% L is set such that the distance between two rays that are 2*pi/L apart
% is one pixel at the outermost radius. Make L even.
% L=ceil(2*pi/atan(2/szvol(1)));
% if mod(L,2)==1 % Make n_theta even
%     L=L+1;
% end
L=360;

% Compute polar Fourier transform of the processed projections
n_r=ceil(szvol(1)/2);
log_message('Start computing polar Fourier transforms of input projections. Using n_r=%d L=%d.',n_r,L);
projs_fname_pft=tempmrcname;
cryo_pft_outofcore(projs_fname_normalized,projs_fname_pft,n_r,L);
log_message('Computing polar Fourier transform done');

projs_fname_pft_n=tempmrcname;
log_message('Start normalizing Fourier transform of input projections (cryo_raynormalize)');
cryo_raynormalize_outofcore(projs_fname_pft,projs_fname_pft_n);
log_message('Normalizing done');
delete(projs_fname_pft);

max_shift=round(szprojs(1)*0.1);
shift_step=0.5;

R1=zeros(3,3,szprojs(3),maxiter);
shift1=zeros(2,szprojs(3),maxiter);
R_refined1=zeros(3,3,szprojs(3),maxiter);
shifts_refined1=zeros(2,szprojs(3),maxiter);
errs1=zeros(szprojs(3),4,maxiter);

iter=1;
roterr=10*tol;
cr=0;

while iter<=maxiter && roterr>tol && cr<cr_threshold
    
    log_message('******** Starting iteration %d/%d ********',iter,maxiter);
    
    log_message('Start orienting projections with respect to reference volume');
    t_orient=tic;
    [R1(:,:,:,iter),shift1(:,:,iter)]=cryo_orient_projections_gpu_outofcore_worker(projs_fname_pft_n,vol,Nrefs,...
        max_shift,shift_step,[],1);
    t_orient=toc(t_orient);
    log_message('Finished orienting projections. Took %7.2f seconds',t_orient);
    
%     log_message('Reconstructing from the projections and their etimated orientation parameters');
%     n=size(projs,1);
%     [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projs,R1,-shift1.', 1e-8, 100, zeros(n,n,n));
%     ii1=norm(imag(v1(:)))/norm(v1(:));
%     log_message('Relative norm of imaginary components = %e\n',ii1);
%     v1=real(v1);
%     
%     fname=sprintf('%s_%d.mrc',map_out_step1,iter);
%     WriteMRC(v1,1,fname);
    
    %poolreopen(8);
    log_message('Start refining orientations');
    [R_refined1(:,:,:,iter),shifts_refined1(:,:,iter),errs1(:,:,iter)]=cryo_refine_orientations_outofcore(...
        projs_fname_pft_n,1,vol,R1(:,:,:,iter),shift1(:,:,iter),1,-1);
    log_message('Finished refining orientations');

    log_message('Start reconstructing from the projections and their refined orientation parameters');
    n=size(projs,1);
    rotations=R_refined1(:,:,:,iter);
    dx=shifts_refined1(:,:,iter);
    [ v1_refined, ~, ~ ,~, ~, ~] = recon3d_firm_outofcore(projs_fname_normalized,rotations,-dx.', 1e-8, 100, zeros(n,n,n));
    log_message('Finished reconstruction');    

    ii1=norm(imag(v1_refined(:)))/norm(v1_refined(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1_refined=real(v1_refined);
    
    fname=sprintf('%s_%d',map_out_step2,iter);
    fname=addExtIfNeeded(fname,'.mrc'); % Add extension if needed.
    WriteMRC(v1_refined,1,fname);
    log_message('Saved %s',fname);
    
    fname=addExtIfNeeded(mat_out,'.mat');
    save(fname)
    log_message('Saved %s',fname);
    
    
    vol=v1_refined;
    
    if iter>=2
        %dd=diag(dist_between_rot(R_refined1(:,:,:,iter),R_refined1(:,:,:,iter-1)))/pi*180;
        dd=rot_dist(R_refined1(:,:,:,iter),R_refined1(:,:,:,iter-1))/pi*180;
        dd=sort(dd);
        log_message('Statistics for iteration %d:',iter)
        log_message('\t Percentiles of assignment errors (in degrees):');
        str=sprintf('%4d\t',10:10:100);
        log_message('\t%s',str);
        str=sprintf('%4.2f\t',dd(floor((10:10:100)/100*numel(dd))));
        log_message('\t %s',str);
        roterr=max(dd);
        
        str=sprintf('%s_%d',map_out_step2,iter-1);
        str=addExtIfNeeded(str,'.mrc');
        vol_prev=ReadMRC(str);
        cr=corr(vol(:),vol_prev(:));
        log_message('\t Correlation of reconstruction of iteration %d with iteration %d is %4.2f',...
            iter,iter-1,cr);
    end
        
    iter=iter+1;
end

% Remove temporary files
delete(projs_fname_normalized);
delete(projs_fname_pft_n)