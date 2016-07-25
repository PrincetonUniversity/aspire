vol1=ReadMRC('/home/yoel/tmp/80S_89/vol_group1.mrc');
projs=ReadMRC('/home/yoel/tmp/80S_89/phaseflipped_downsampled_prewhitened_group1.mrc');
[Rests,dxs]=cryo_orient_projections_gpu(projs,vol1,[],[],0,0);

n=size(projs,1);
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( projs,Rests,dxs.', 1e-6, 30, zeros(n,n,n));
ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v1=real(v1);
WriteMRC(v1,1,'/home/yoel/tmp/80S_89/v1.mrc');

vol2=ReadMRC('/home/yoel/tmp/80S_89/vol_group2_aligned.mrc');
projs=ReadMRC('/home/yoel/tmp/80S_89/phaseflipped_downsampled_prewhitened_group2.mrc');
[Rests,dxs]=cryo_orient_projections_gpu(projs,vol2,[],[],0,0);

n=size(projs,1);
[ v2, ~, ~ ,~, ~, ~] = recon3d_firm( projs,Rests,dxs.', 1e-6, 30, zeros(n,n,n));
ii1=norm(imag(v2(:)))/norm(v2(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
v2=real(v2);
WriteMRC(v2,1,'/home/yoel/tmp/80S_89/v2.mrc');

plotFSC(vol1,vol2,0.143,5.77);
plotFSC(v1,v2,0.143,5.77);

save tmp
