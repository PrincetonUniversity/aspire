% A call to impl with only mandatory variables
cache_file_name = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/cn/cn_cache_points1000_ntheta360_res1.mat';
verbose = 0;
%%
% n_symm = 5;
% empiar_code_string = 10089;
% % mrc_stack_file = '/scratch/yoel/denoised/10089/denoised_group1.mrcs';


% %%
% n_symm = 4;
% empiar_code_string = '10081';
% % emdb_code = 8511;
% % produces nice results but resolution against the emdb (0.5 criterion) is
% % 35A. Need to try older C4
% mrc_stack_file = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/cn/datasets_c4/10081/averages_nn50_group1.mrc';
% % down_siz = 89; 
% % pixA = 1.3*256/down_siz;

% % %
% n_symm = 11;
% empiar_code_string = '10063';
% emdb_code = 6458;
% mrc_stack_file = '/home/yoel/scratch/10063/aspire/ring11/averages_nn50_group1.mrc';
% down_siz = 89;
% pixA = 0.86*448/down_siz;




% %
% %     instack = '10004_averages_nn100_selected_small.mrcs';
% n_symm = 3;
% empiar_code_string = '10004';
% emdb_code = 2484;
% mrc_stack_file = '/home/yoel/scratch/10004/averages_nn100_selected.mrcs';
% mask_radius_perc = 50;
% down_siz = 89;
% pixA = 2.16*127/down_siz;
% %     xml_file = '/home/yoel/scratch/10004/aspire/EMPIAR1004symC3.xml';


n_symm = 7;
empiar_code_string = 'GroEL';
mrc_stack_file = '/scratch/yoel/ranz/aspire89/averages_nn100_group2.mrc';

recon_folder = '/home/gabip/matlabProjects/aspire/aspire/development/abinitio/cn/results';
recon_mrc_fname = fullfile(recon_folder,sprintf('%s_out.mrc',empiar_code_string));
recon_mat_fname = fullfile(recon_folder,sprintf('%s_out.mat',empiar_code_string));

log_fname = fullfile(recon_folder,'log.txt');
open_log(log_fname);

cryo_abinitio_Cn(cache_file_name,n_symm,mrc_stack_file,recon_mrc_fname,recon_mat_fname,verbose);

close_log();


% cryo_abinitio_C4('/home/yoel/scratch/fred/aspire89/averages_nn100_group1.mrc','fred_c4_89_nn100_grp1_2000_linspace.mrc','fred_c4_89_nn100_grp1_2000_linspace.mat',2000);
% clear; close all;
% cryo_abinitio_C4('/home/yoel/scratch/fred/aspire89/averages_nn100_group2.mrc','fred_c4_89_nn100_grp2_2000_linspace.mrc','fred_c4_89_nn100_grp2_2000_linspace.mat',2000);
% clear; close all;
% 
% vol1 = ReadMRC('fred_c4_89_nn100_grp1_2000_linspace.mrc');
% vol2 = ReadMRC('fred_c4_89_nn100_grp2_2000_linspace.mrc');
% pixA = 3.3;
% [Rest,estdx,vol2aligned] = cryo_align_densities_C4(vol1,vol2,pixA,1,[],0,50);
% [resA,h] = plotFSC(vol1,vol2aligned,0.143,pixA);
% 



%%
% n_symm = 4;
% mrc_stack_file = '/scratch/yoel/fred/aspire89/averages_nn50_group1.mrc';
% empiar_code_string = 'fred';
% down_siz = 89;
% pixA = 192*1.533/down_siz;


% vol_1_mat = cryo_fetch_emdID(emdb_code);
% vol_1 = ReadMRC(vol_1_mat);
% vol_1 = cryo_downsample(vol_1,down_siz,0);
% 
% vol_2 = ReadMRC(recon_mrc_fname);
% vol_2 = cryo_downsample(vol_2,down_siz,0);
% 
% [~,~,vol2aligned] = cryo_align_densities_C4(vol_1,vol_2,pixA,1);
% 
% [resA,h] = plotFSC(vol_1,vol2aligned,0.5,pixA);