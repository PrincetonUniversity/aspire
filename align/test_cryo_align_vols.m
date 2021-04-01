%% test_cryo_align_vols
% test for cryo_align_vols_rot_symmetry

%% example of cyclic molecule:
mapfile=cryo_fetch_emdID(10280);   % C1 symmetry
%mapfile=cryo_fetch_emdID(0659);   % C1 symmetry
vol=ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
symmetry = 'C';
n_symmetry = 1;
G = cryo_cyclic_symmetry_group(n_symmetry);


%% example of C6 molecule:
%mapfile=cryo_fetch_emdID(10477);   % C6 symmetry
mapfile=cryo_fetch_emdID(0825);

vol=ReadMRC(mapfile);
sz_vol = 127;
vol = cryo_downsample(vol,sz_vol,0);
symmetry = 'C';
n_symmetry = 6;
G = cryo_cyclic_symmetry_group(n_symmetry);
%G = findSymElements(vol,6);

%% example of Dihedral molecule:
symmetry = 'D';
n_symmetry = 3;
mapfile = cryo_fetch_emdID(20656);
%mapfile = cryo_fetch_emdID(9203); 
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
G = cryo_dihedral_symmetry_group2(3);

%% example of cubic symmetry:
symmetry = 'T';
n_symmetry = 1;
mapfile = cryo_fetch_emdID(4179); 
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
G = cryo_T_symmetry_group();

%% example of cubic symmetry:
symmetry = 'I';
n_symmetry = 1;
%mapfile = cryo_fetch_emdID(4550);
mapfile = cryo_fetch_emdID(0069);
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
G = cryo_I_symmetry_group();

%% example of cubic symmetry:
symmetry = 'O';
n_symmetry = 1;
%mapfile = cryo_fetch_emdID(22346);
mapfile = cryo_fetch_emdID(11121);
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
G = cryo_O_symmetry_group();

%% D2-beta gal
symmetry = 'D';
n_symmetry = 2;
mapfile = cryo_fetch_emdID(7770);
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
n = 4;
G = cryo_dihedral_symmetry_group(2);

%%
figure
view3d(vol)
title('tested reference volume')



%R = rand_rots(1);
%%% rotate vol:
%vol_ref = fastrotate3d(vol,R);
%figure
%view3d(vol_ref)
%title('reference volume')

true_R = rand_rots(1);
vol_rotated = fastrotate3d(vol,true_R);


%%% add reflection:
vol_rotated = flip(vol_rotated,3);


%%% shift vol:
vol_rotated = reshift_vol(vol_rotated,[-5 0 0]);
figure
view3d(vol_rotated)
title('rotated and shifted volume')


%% Add noise to volume:
noise = 1;
SNR = 0.1;

sigma = sqrt(1/SNR);
%n1 = sigma*randn(size(vol))./sqrt(numel(vol(:)));
n2 = sigma*randn(size(vol_rotated))./sqrt(numel(vol_rotated(:)));
%vol_Noisy = vol+n1;
vol_rotated = vol_rotated + n2;

figure
view3d(vol_rotated)
title('rotated and shifted volume + noise')

%% align volumes:
N_projs = 30;      % for the alignment.

t1 = tic;
[R_est,estdx,vol_aligned,corr_R] = cryo_align_vols(symmetry,n_symmetry,vol,vol_rotated,[],true_R,G);
T_align = toc(t1);
log_message('The alignment took %5.2f seconds.',T_align);






