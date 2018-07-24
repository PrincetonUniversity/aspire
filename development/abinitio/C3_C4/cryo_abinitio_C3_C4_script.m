% instack = '/scratch/yoel/denoised/10081/denoised_group1.mrcs';

%%
% %     instack = '10004_averages_nn100_selected_small.mrcs';
% n_symm = 3;
% instack = '/home/yoel/scratch/10004/averages_nn100_selected.mrcs';
% mask_radius_perc = 50;
% recon_folder = './results_tmp_C3';
% empiar_code_string = '10004';
% emdb_code = 2484;
% down_siz = 89;
% pixA = 2.16*127/down_siz;
% %     xml_file = '/home/yoel/scratch/10004/aspire/EMPIAR1004symC3.xml';
% 
% %%
% n_symm = 4;
% instack = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C3_C4/datasets_C4/10081/averages_nn50_group1.mrc';
% %     instack = '10081_denoised_group1_small.mrcs';
% mask_radius_perc = 30;
% recon_folder = './results_tmp_C4';
% empiar_code_string = '10081';
% emdb_code = 8511;
% down_siz = 89;
% pixA = 1.3*255/down_siz;

%%
n_symm = 4;
instack = '/scratch/yoel/fred/aspire89/averages_nn50_group1.mrc';
recon_folder = './results_tmp_C4';
mask_radius_perc = 40;
empiar_code_string = 'fred';
down_siz = 89;
pixA = 192*1.533/down_siz;

%%
outvol = fullfile(recon_folder,sprintf('%s_out.mrc',empiar_code_string));

%% option 1: using minimal set of input variables
% cryo_abinitio_C3_C4(n_symm,instack,outvol);

%% option 2: using all input variables
outmat = fullfile(recon_folder,sprintf('%s_out.mat',empiar_code_string));
max_shift_perc = 15;
shift_step = 0.5;
n_r_perc = 50;

n_theta = 360;

log_fname = fullfile(recon_folder,'log.txt');
open_log(log_fname);
cryo_abinitio_C3_C4(n_symm,instack,outvol,outmat,max_shift_perc,shift_step,n_r_perc,mask_radius_perc,n_theta);


vol_1_mat = cryo_fetch_emdID(emdb_code);
vol_1 = ReadMRC(vol_1_mat);
vol_1 = cryo_downsample(vol_1,down_siz,0);

vol_2 = ReadMRC(outvol);
% vol_2 = ReadMRC('./results_C4/vol_10081_nn50_grp1_100.mrc');
% vol_2 = ReadMRC('./results_C3/c3_10004_nn100_selected_ims3500.mrc');
vol_2 = cryo_downsample(vol_2,down_siz,0);

[~,~,vol2aligned] = cryo_align_densities_C4(vol_1,vol_2,pixA,1,[],0,100);

[resA,h] = plotFSC(vol_1,vol2aligned,0.5,pixA);


close_log();