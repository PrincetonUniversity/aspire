
vol_init_mat = '/tmp/tpac3abed0_342e_41ac_89e5_d2e8add6e08f/pub/databases/emdb/structures/EMD-2660/map/emd_2660.map';
vol_init = ReadMRC(vol_init_mat);
vol_init = cryo_downsample(vol_init,89);
clear vol_init_mat

vol=genD2fromVol(vol_init);
projs = cryo_project(vol,permute(test_grid,[2,1,3]));
projs = permute(projs,[2,1,3]);
viewstack(projs(:,:,1:2),2,1,0);
