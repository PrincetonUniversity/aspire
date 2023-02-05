% Test the function cryo_abinitio_TO using simulated noisy data with T symmetry
% 
% Generate output files for T symmetry.
% The figures of the paper are generated from these output files using the function plot_results_T.m
%
% Yoel Shkolnisky, October 2022.

out_dir = './results_T';

mapfile = cryo_fetch_emdID(10835); % Download 3D map file with O symmetry
volref=ReadMRC(mapfile);           % Load the map file
sz_orig=size(volref,1);
pixA=0.639;
sz_ds=129;
volref_ds=cryo_downsample(volref,sz_ds);  % Downsample map   
WriteMRC(volref_ds(1:end-1,1:end-1,1:end-1),pixA*sz_orig/sz_ds,...
    fullfile(out_dir,'volref.mrc'));


for Nprojs=[25, 50, 100, 200]
    initstate;
    rots=rand_rots(Nprojs); % Generate random rotations
    projs=cryo_project(volref_ds,rots);
    projs=permute(projs,[2 1 3]);
    
    snr_list = [1000, 1, 1/2, 1/4, 1/8];
    
    for k=1:numel(snr_list)
        snr=snr_list(k);
        noisy_projs = cryo_addnoise(projs,snr,'gaussian');
        projs_fname=sprintf('noisy_projs_%d_%d.mrcs',Nprojs,k);
        projs_fname=fullfile(out_dir,projs_fname);
        
        vol_fname=sprintf('vol_%d_%d.mrc',Nprojs,k);
        vol_fname=fullfile(out_dir,vol_fname);
        
        params_fname=sprintf('params_%d_%d.mat',Nprojs,k);
        params_fname=fullfile(out_dir,params_fname);
        
        WriteMRC(noisy_projs,pixA*sz_orig/sz_ds,projs_fname); % Save simulated projections
        
        t_start = tic;
        cryo_abinitio_TO('T',projs_fname,vol_fname,...
            'T_symmetry_2084_candidates_cache_080222.mat',params_fname);
        execution_time = toc(t_start);
        timing_fname=sprintf('timing_%d_%d.txt',Nprojs,k);
        fid = fopen(fullfile(out_dir,timing_fname),'w');
        fprintf(fid,'%d\n',execution_time);
        fclose(fid);

        % Check reconstruction
        figure;
        vol2 = ReadMRC(vol_fname);
        [bestR,bestdx,~,vol2aligned,bestcorr]=cryo_align_vols(volref_ds(1:end-1,1:end-1,1:end-1),vol2,1);
        
        % Save aligned volume
        vol_fname=sprintf('vol_%d_%d_aligned.mrc',Nprojs,k);
        vol_fname=fullfile(out_dir,vol_fname);
        WriteMRC(vol2aligned,pixA*sz_orig/sz_ds,vol_fname);
        
        plotFSC(volref_ds(1:end-1,1:end-1,1:end-1),vol2aligned,0.5,pixA*sz_orig/sz_ds,1);
        savefig(fullfile(out_dir,sprintf('fsc_%d_%d.fig',Nprojs,k)));
    end
    
end
