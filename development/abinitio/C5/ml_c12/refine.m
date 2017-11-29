% 
% pixA = 0.86;
% 
% vol1 = ReadMRC('/scratch/gabip/10063/gabi_abinitio/ring11/ring11/job900_it025_class001.mrc');
% vol2 = ReadMRC('/scratch/gabip/10063/emdb6458/_it025_class001.mrc');
% 
% [~,~,vol2aligned] = cryo_align_densities(vol1,vol2,pixA,1);
% 
% [resA,h] = plotFSC(vol1,vol2aligned,0.143,pixA);
% 
% 
% parpool(2);
% down_siz = 65;
% pixA = 0.86*448/down_siz;
% 
% vol_1_mat = cryo_fetch_emdID(6458);
% vol_1 = ReadMRC(vol_1_mat);
% vol_1 = cryo_downsample(vol_1,down_siz,0);
% 
% vol_2 = ReadMRC('/scratch/gabip/10063/gabi_abinitio/10063_448pix_c12.mrc');
% vol_2 = cryo_downsample(vol_2,down_siz,0);
% 
% [~,~,vol2aligned] = cryo_align_densities_C4(vol_1,vol_2,pixA,1);
% 
% [resA,h] = plotFSC(vol_1,vol2aligned,0.5,pixA);

parpool(2);
down_siz = 65;
pixA = 0.86*448/down_siz;

vol_1 = ReadMRC('/scratch/gabip/10063/emdb6458/emdb6458_abinitio_lowpassed.mrc');
vol_1 = cryo_downsample(vol_1,down_siz,0);

vol_2 = ReadMRC('/scratch/gabip/10063/emdb6458/refined2/_it025_class001.mrc');
vol_2 = cryo_downsample(vol_2,down_siz,0);

[~,~,vol2aligned] = cryo_align_densities_C4(vol_1,vol_2,pixA,1);

[resA,h] = plotFSC(vol_1,vol2aligned,0.5,pixA);



% vol_1_mat = cryo_fetch_emdID(6458);
% vol_1 = ReadMRC(vol_1_mat);
% % code for low_passed %%
% vol_1_lowpassed = GaussFilt(vol_1,0.2);
% WriteMRC(vol_1_lowpassed,1,'/scratch/gabip/10063/emdb6458/emdb6458_abinitio_lowpassed.mrc');



% [~,~,vol2aligned] = cryo_align_densities_C4(vol_1,vol_2,pixA,1);
% 
% [resA,h] = plotFSC(vol_1,vol2aligned,0.5,pixA);



