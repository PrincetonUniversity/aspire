
OUTPUTDIR='/tmp/clml';

for N=[500,1000]
    % Read projections
    projs=ReadMRC('~/tmp/80S_89/averages_nn50_group1.mrc');
    projs=projs(:,:,1:N);

    % Mask projections
    K=size(projs,3);
    mask_radius=round(size(projs,1)*0.45);
    [masked_projs,~]=mask_fuzzy(projs,mask_radius);

    % Estimate noise power spectrum
    spnoise_pft= cryo_noise_estimation_pfs(projs,n_r,n_theta);
    spnoise_pft(spnoise_pft<max(spnoise_pft/10))=1;
    noise_cov=diag(spnoise_pft);

    % Polar Fourier transform
    n_theta=360;
    n_r=ceil(size(projs,1)*0.5);
    [pf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections

    % Find common lines from projections using maximum
    % likelihood
    pf=cryo_raynormalize(reshape(pf,n_r,n_theta*N));
    pf=reshape(pf,n_r,n_theta,N);
    
    for max_shift=[10, 15]
        open_log(0);        
        shift_step=1;
        
        for max_iterations=[4,6,8]
            M=26000;
            KNN=50;
            [clstack,~,shift_equations,shift_equations_map,~]=...
                cryo_clmatrix_ml_gpu(pf,M,KNN,noise_cov,max_iterations,...
                1,max_shift,shift_step);
            
            S=cryo_syncmatrix_vote(clstack,n_theta);
            rotations=cryo_syncrotations(S);
            
            clerr=syncconsistency(rotations,clstack,n_theta);
            figure(1);clf;hist(clerr(:,3),360);
            figname=sprintf('clerr_N%d_ms%d_mi%d.eps',N,max_shift,max_iterations);
            print(fullfile(OUTPUTDIR,figname),'-depsc');
            
            reloadname=sprintf('reload_real_data_eee_N%d_ms%d_mi%d',...
                N,max_shift,max_iterations);
            save(fullfile(OUTPUTDIR,reloadname));
        end
    end
end