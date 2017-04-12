function cryo_refine_map_iter(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out,maxiter,tol)
%
% XXX Add convergence
% XXX Add CTF
% XXX Do we need filtering during assignment of orientations?
% XXX Create function cryo_mask(vol,radius,isstack)

if ~exist('maxiter','var')
    maxiter=5;
end

if ~exist('tol','var')
    tol=10;
end

log_message('Using maxiter=%d tol=%4.2f degrees',maxiter,tol);
log_message('Iterations stop when maxiter reached or the rotations of the');
log_message('current iteration deviate from the rotations of the previous');
log_message('iteatrion by a maximum of tol degrees');
    

vol=ReadMRC(vol_mrc);
projs=ReadMRC(projs_mrc);

szvol=size(vol);
szprojs=size(projs);

projs=cryo_globalphaseflip(projs);
vol=cryo_globalphaseflip_vol(vol);


if szvol(1)~=szprojs(1)
    vol=cryo_downsample(vol,[szprojs(1),szprojs(1),szprojs(1)]);
end

R1=zeros(3,3,szprojs(3),maxiter);
shift1=zeros(2,szprojs(3),maxiter);
R_refined1=zeros(3,3,szprojs(3),maxiter);
shifts_refined1=zeros(2,szprojs(3),maxiter);
errs1=zeros(szprojs(3),4,maxiter);

iter=1;
roterr=10*tol;
cr=0;

while iter<=maxiter && roterr>tol && cr<0.8
    
    log_message('******** Starting iteration %d/%d ********',iter,maxiter);
    
    t_orient=tic;
    %[R1,shift1]=cryo_orient_projections_gpu_2(projs,vol,-1,[],1,1,4);
    [R1(:,:,:,iter),shift1(:,:,iter)]=cryo_orient_projections_gpu(projs,vol,-1,[],1,1);
    %[R1,shift1]=cryo_orient_projections(projs,vol,[],[],1,1);
    t_orient=toc(t_orient);
    
    log_message('First step of refinement took %7.2f seconds',t_orient);
    
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
    [R_refined1(:,:,:,iter),shifts_refined1(:,:,iter),errs1(:,:,iter)]=cryo_refine_orientations(...
        projs,vol,R1(:,:,:,iter),shift1(:,:,iter),1,-1);
    
    log_message('Reconstructing from the projections and their refined orientation parameters');
    n=size(projs,1);
    rotations=R_refined1(:,:,:,iter);
    dx=shifts_refined1(:,:,iter);
    [ v1_refined, ~, ~ ,~, ~, ~] = recon3d_firm( projs,rotations,-dx.', 1e-8, 100, zeros(n,n,n));
    ii1=norm(imag(v1_refined(:)))/norm(v1_refined(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1_refined=real(v1_refined);
    
    fname=sprintf('%s_%d.mrc',map_out_step2,iter);
    WriteMRC(v1_refined,1,fname);

    fname=sprintf('%s.mat',mat_out);
    save(fname)
    
    vol=v1_refined;

    
    if iter>=2
        dd=diag(dist_between_rot(R_refined1(:,:,:,iter),R_refined1(:,:,:,iter-1)))/pi*180;
        dd=sort(dd);
        log_message('Percentiles of assignment errors (in degrees):');
        str=sprintf('%4d\t',10:10:100);
        log_message('%s',str);
        str=sprintf('%4.2f\t',dd(floor((10:10:100)/100*numel(dd))));
        log_message('%s',str);
        roterr=max(dd);
        
        str=sprintf('%s_%d.mrc',map_out_step2,iter-1);
        vol_prev=ReadMRC(str);
        cr=corr(vol(:),vol_prev(:));
        log_message('Correlation of reconstruction of iteration %d with iteration %d is %4.2f',...
            iter,iter-1,cr);
    end
        
    iter=iter+1;
end

