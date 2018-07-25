% CTFdata = readSTAR('/home/yoel/scratch/10063/particles/ring10.star');
% CTFdata = addfieldtoSTARdata(CTFdata,'pixA',1.72);
% writeSTAR(CTFdata,'ring10_new.star');

vol = ReadMRC('10063_nn50_ring11_grp1_ims1000_ml_reinfC11_filt.mrc');
vol_out = cryo_downsample(vol,448);
WriteMRC(vol_out,1,'10063_nn50_ring11_grp1_ims1000_ml_reinfC11_filt_upsmpld.mrc');
