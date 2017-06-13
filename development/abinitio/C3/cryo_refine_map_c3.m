function cryo_refine_map_c3(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out,maxiter,tol)

if ~exist('maxiter','var')
    maxiter=5;
end

if ~exist('tol','var')
    tol=10;
end
cr_threshold=0.8;

log_message('Using maxiter=%d, tol=%4.2f degrees, correlation threshold=%d;',maxiter,tol,cr_threshold);
log_message('Iterations stop when maxiter reached or the rotations of the current iteration deviate from the rotations of the previous');
log_message('iteatrion by a less than tol degrees, or the correlation between volumes of two consecutive iterations is at least as specified.');
    
log_message('Loading initial volume %s',vol_mrc');
vol=ReadMRC(vol_mrc);
log_message('Volume loaded');
log_message('vol MD5 %s',MD5var(vol));

log_message('Loading projections %s',projs_mrc);
projs=ReadMRC(projs_mrc);


log_message('***************************************');
log_message('***************************************');
log_message('***************************************');
log_message('REMOVE THIS LINE');
log_message('SAMPLING 2000 IMAGES');
projs = projs(:,:,randperm(size(projs,3),2000));
log_message('***************************************');
log_message('***************************************');
log_message('***************************************');
log_message('REMOVE THIS LINE');


log_message('Projections loaded');
log_message('projs MD5 %s',MD5var(projs));

szvol=size(vol);
szprojs=size(projs);

log_message('Applying global phase-flip to initial volume');
vol=cryo_globalphaseflip_vol(vol);
log_message('Volume phase-flipped');
log_message('vol MD5 %s',MD5var(vol));

log_message('Applying global phase-flip projections');
projs=cryo_globalphaseflip(projs);
log_message('Projections phase-flipped');
log_message('projs MD5 %s',MD5var(projs));

log_message('Dimensions of volume %dx%dx%d]',szvol(1),szvol(2),szvol(3));
log_message('Dimensions of projections %dx%d (%d projections)',szvol(1),szvol(2),szvol(3));
if szvol(1)~=szprojs(1)
    log_message('Volume and projections have differnt dimensions. Resampling initial volume to size %dx%dx%d',...
        szprojs(1),szprojs(1),szprojs(1));
    vol=cryo_downsample(vol,[szprojs(1),szprojs(1),szprojs(1)]);
    log_message('Initial volume resampled');
end
log_message('vol MD5 %s',MD5var(vol));

szvol=size(vol); % Dimensions of volume may have changed due to resampling
% Preprocess projections and referece volume.
log_message('Start preprocessing projections');
n=szvol(1);
projs=cryo_mask(projs,1,floor(0.45*n),floor(0.05*n)); % Mask projections
log_message('Preprocessing done');
log_message('projs MD5 %s',MD5var(projs));

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
projs_hat=cryo_pft(projs,n_r,L,'single');
projs_hat=single(projs_hat);
log_message('Computing polar Fourier transform done');
log_message('projs_hat MD5 %s',MD5var(projs_hat));


log_message('Start normalizing Fourier transform of input projections (cryo_raynormalize)');
projs_hat=cryo_raynormalize(projs_hat);
log_message('Normalizing done');
log_message('projs_hat MD5 %s',MD5var(projs_hat));

max_shift=round(size(projs,1)*0.1);
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
    [R1(:,:,:,iter),shift1(:,:,iter)]=cryo_orient_projections_gpu_worker(projs_hat,vol,Nrefs,...
        max_shift,shift_step,[],1);
    t_orient=toc(t_orient);
    log_message('Finished orienting projections. Took %7.2f seconds',t_orient);
    log_message('R1 MD5 %s',MD5var(R1));
    log_message('shift1 MD5 %s',MD5var(shift1));

    
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
    log_message('projs MD5 %s',MD5var(projs));
    log_message('vol MD5 %s',MD5var(vol));
    log_message('projs_hat MD5 %s', MD5var(projs_hat));
    log_message('Start refining orientations');
    [R_refined1(:,:,:,iter),shifts_refined1(:,:,iter),errs1(:,:,iter)]=cryo_refine_orientations(...
        projs_hat,1,vol,R1(:,:,:,iter),shift1(:,:,iter),1,-1);
    log_message('Finished refining orientations');
    log_message('R_refined1 MD5 %s',MD5var(R_refined1));
    log_message('shifts_refined1 MD5 %s',MD5var(shifts_refined1));
    log_message('errs1 MD5 %s',MD5var(errs1));

    log_message('Start reconstructing from the projections and their refined orientation parameters');
    n=size(projs,1);
    rotations=R_refined1(:,:,:,iter);
    dx=shifts_refined1(:,:,iter);
    log_message('projs MD5 %s',MD5var(projs));
    log_message('rotations MD5 %s',MD5var(rotations));
    log_message('dx MD5 %s',MD5var(dx));
    v1_refined = reconstruct(projs,rotations,n_r,L,max_shift,shift_step);
%     [ v1_refined, ~, ~ ,~, ~, ~] = recon3d_firm( projs,rotations,-dx.', 1e-8, 100, zeros(n,n,n));
    log_message('Finished reconstruction');
    log_message('v1_refined MD5 %s',MD5var(v1_refined));

    ii1=norm(imag(v1_refined(:)))/norm(v1_refined(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1_refined=real(v1_refined);
    
    fname=sprintf('%s_%d.mrc',map_out_step2,iter);
    WriteMRC(v1_refined,1,fname);

    fname=sprintf('%s.mat',mat_out);
    save(fname)
    
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
        
        str=sprintf('%s_%d.mrc',map_out_step2,iter-1);
        vol_prev=ReadMRC(str);
        cr=corr(vol(:),vol_prev(:));
        log_message('\t Correlation of reconstruction of iteration %d with iteration %d is %4.2f',...
            iter,iter-1,cr);
    end
        
    iter=iter+1;
end
