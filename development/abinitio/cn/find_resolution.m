% 1.04 is from empiar. 256 is the original pixel size (from xml), 129 is the downsamples size (from xml), but in code it was downsampled again to 89
im_sz = 89;
pix_size = 2.7819; 
mapFile = cryo_fetch_emdID(2981);
vol_ref = ReadMRC(mapFile);
vol_ref_downsmpl = cryo_downsample(vol_ref,[im_sz,im_sz,im_sz]);

vol2 = ReadMRC('c3_10097_nn50_cnstrst_grp1_ims1000.mrc'); 
% vol2 = ReadMRC('10038_refspick_nn10_grp1_ims1000.mrc');
[~,~,vol2_aligned] = cryo_align_densities_C2(vol_ref_downsmpl,vol2,pix_size,1,[],0,100);
[resA, h] = plotFSC(vol_ref_downsmpl,vol2_aligned,0.143,pix_size);