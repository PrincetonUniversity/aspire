mapFile = cryo_fetch_emdID(6267);
vol_ref = ReadMRC(mapFile);
vol_ref_down_sam = cryo_downsample(vol_ref,[129,129,129]);

vol2 = ReadMRC('vol_10024_grp1_3000.mrc');
[~,~,vol2aligned] = cryo_align_densities_C4(vol_ref_down_sam,vol2,2.82,1,[],0,100);
[resA,h] = plotFSC(vol_ref_down_sam,vol2aligned,0.143,2.82);



vol1=ReadMRC('results/10081/vol_10081_nn50_grp1_5000.mrc');
vol2=ReadMRC('results/10081/vol_10081_nn50_grp2_5000.mrc');

[~,~,vol2aligned] = cryo_align_densities_C4(vol1,vol2,2.5698,1,[],0,100);
[resA,h] = plotFSC(vol1,vol2aligned,0.143,2.5698);