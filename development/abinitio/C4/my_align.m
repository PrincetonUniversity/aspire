vol1 = ReadMRC('vol_10081_nn50_grp1_5000.mrc');
vol2 = ReadMRC('vol_10081_nn50_grp2_5000.mrc');
[Rest,estdx,vol2aaligned]=cryo_align_densities_C4(vol1,vol2,2.5698,1,[],0,50);
[resA,h]=plotFSC(vol1,vol2aligned,0.143,2.5698);

