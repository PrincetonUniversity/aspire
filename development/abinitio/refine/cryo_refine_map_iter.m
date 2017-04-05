function cryo_refine_map_iter(projs_mrc,vol_mrc,map_out_step1,map_out_step2,mat_out)

vol=ReadMRC(vol_mrc);
projs=ReadMRC(projs_mrc);
Niter=3;

szvol=size(vol);
szprojs=size(projs);

projs=cryo_globalphaseflip(projs);
vol=cryo_globalphaseflip_vol(vol);


if szvol(1)~=szprojs(1)
    vol=cryo_downsample(vol,[szprojs(1),szprojs(1),szprojs(1)]);
end

for iter=1:Niter
    log_message('******** Starting iteration %d/%d ********',iter,Niter);
    
    t_orient=tic;
    %[R1,shift1]=cryo_orient_projections_gpu_2(projs,vol,-1,[],1,1,4);
    [R1,shift1]=cryo_orient_projections_gpu(projs,vol,-1,[],1,1);
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
    [R_refined1,shifts_refined1,errs1]=cryo_refine_orientations(...
        projs,vol,R1,shift1,1,-1);
    
    log_message('Reconstructing from the projections and their refined orientation parameters');
    n=size(projs,1);
    [ v1_refined, ~, ~ ,~, ~, ~] = recon3d_firm( projs,R_refined1,-shifts_refined1.', 1e-8, 100, zeros(n,n,n));
    ii1=norm(imag(v1_refined(:)))/norm(v1_refined(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1_refined=real(v1_refined);
    
    fname=sprintf('%s_%d.mrc',map_out_step2,iter);
    WriteMRC(v1_refined,1,fname);

    fname=sprintf('%s_%d.mat',mat_out,iter);
    save(fname)
    
    vol=v1_refined;
end

