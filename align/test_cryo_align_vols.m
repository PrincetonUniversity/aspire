%% test_cryo_align_vols


%% example of cyclic molecule:
%mapfile = cryo_fetch_emdID(10280);   % C1 symmetry
%sym = 'C1';
%% example of C6 molecule:
%mapfile=cryo_fetch_emdID(10477);   % C6 symmetry
%mapfile = cryo_fetch_emdID(0825);
%sym = 'C6';
%%
%mapfile=cryo_fetch_emdID(11516);
%sym = 'C7';
%%
%mapfile=cryo_fetch_emdID(24494);
%mapfile=cryo_fetch_emdID(3528);
%mapfile=cryo_fetch_emdID(22854);
%sym = 'I';
%% D2-beta gal
%mapfile = cryo_fetch_emdID(7770);
%sym = 'D2';
%% 
mapfile = cryo_fetch_emdID(9203); 
sym = 'D3';
%% 
%mapfile = cryo_fetch_emdID(4179);  
%sym = 'T';
%%
%mapfile = cryo_fetch_emdID(22658);
%sym = 'O';
%% from yoel:
%vol1 = ReadMRC('abinitio_10272_ref.mrc');
%vol2 = ReadMRC('abinitio_10272.mrc');
%sym = 'O';
%opt.downsample = 48;
%opt.sym = sym;
%[R_est_yael,estdx_yael,reflect_yael,vol_aligned_yael,corr_R_yael] = cryo_align_vols(vol1,vol2,1,opt);

%%
initstate
%rand('twister', 1337);

%%% Rotate the volume to generate reference volume:
vol = ReadMRC(mapfile);
sz_vol = 129;
vol = cryo_downsample(vol,sz_vol,0);
%[R,~,~] = svd(rand(3)); % Generate random rotation
%vol = fastrotate3d(vol,R);

%figure
%view3d(vol)
%title('tested reference volume')

%%% Rotate the reference volume:
[true_R,~,~] = svd(rand(3)); % Generate random rotation
vol_rotated = fastrotate3d(vol,true_R);

%%% Add reflection:
%vol_rotated = flip(vol_rotated,3);

%%% Reshift the volume:
vol_rotated = reshift_vol(vol_rotated,[-5 0 0]);
%figure
%view3d(vol_rotated)
%title('rotated and shifted volume')
%% Add noise to volume:
%{
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
%}
%% Align volumes:
G = genSymGroup(sym);
n = size(G,3);
%for i = 1:n
%    G(:,:,i) = R*G(:,:,i)*R.';
%end
opt.N_projs = 30;     
opt.G = G;
opt.true_R = true_R;
opt.dofscplot = 1;
opt.downsample = 64;
%opt.sym = sym;

t1 = tic;
[R_est,estdx,reflect,vol_aligned,corr_R] = cryo_align_vols(vol,vol_rotated,1,opt);
T_align = toc(t1);
log_message('The alignment took %5.2f seconds.',T_align);

figure
view3d(vol_aligned)
title('aligned volume')



