vol1=ReadMRC('vol_nn50_nm5000_group1.mrc');
projs1=ReadMRC('phaseflipped_cropped_downsampled_prewhitened_group1.mrc');
t_orient=tic;
[R1,shift1]=cryo_orient_projections_gpu(projs1,vol1,-1,[],1,1);
t_orient=toc(t_orient);

log_message('Reconstructing from the projections and their etimated orientation parameters');
n=size(projs1,1);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projs1,R1,-shift1.', 1e-8, 100, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v1=real(v1);

[R_refined1,shifts_refined1,errs1]=cryo_refine_orientations(...
    projs1,vol1,R1,shift1,1,-1);

log_message('Reconstructing from the projections and their refined orientation parameters');
n=size(projs1,1);
[ v1_refined, ~, ~ ,~, ~, ~] = recon3d_firm( projs1,R_refined1,-shifts_refined1.', 1e-8, 100, zeros(n,n,n));
ii1=norm(imag(v1_refined(:)))/norm(v1_refined(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v1_refined=real(v1_refined);

save tmp1

clear;

vol2=ReadMRC('vol_nn50_nm5000_group2.mrc');
projs2=ReadMRC('phaseflipped_cropped_downsampled_prewhitened_group2.mrc');
t_orient=tic;
[R2,shift2]=cryo_orient_projections_gpu(projs2,vol2,-1,[],1,1);
t_orient=toc(t_orient);

log_message('Reconstructing from the projections and their etimated orientation parameters');
n=size(projs2,1);
[ v2, ~, ~ ,~, ~, ~] = recon3d_firm( projs2,R2,-shift2.', 1e-8, 100, zeros(n,n,n));
ii2=norm(imag(v2(:)))/norm(v2(:));
log_message('Relative norm of imaginary components = %e\n',ii2);
v2=real(v2);

[R_refined2,shifts_refined2,errs2]=cryo_refine_orientations(...
    projs2,vol2,R2,shift2,1,-1);

log_message('Reconstructing from the projections and their refined orientation parameters');
n=size(projs2,1);
[ v2_refined, ~, ~ ,~, ~, ~] = recon3d_firm( projs2,R_refined2,-shifts_refined2.', 1e-8, 100, zeros(n,n,n));
ii2=norm(imag(v2_refined(:)))/norm(v2_refined(:));
log_message('Relative norm of imaginary components = %e\n',ii2);
v2_refined=real(v2_refined);

save tmp2

