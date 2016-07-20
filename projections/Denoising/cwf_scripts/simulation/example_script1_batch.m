% Denoise projections of 80s ribosome (with white noise and CTF added).
[fpath_user]=fmtinput('Please enter the destination path for denoised images. The default destination is your home directory','~/.','%s');
vol=load('cleanrib.mat');

K=1000;
[projections] = cryo_gen_projections(35,K,1,0);

nbatch=10;

N_images=K;
g_projections=projections(:,:,1:N_images);
g_projections=mask_corners(g_projections);
clear projections


hatI=cfft2(g_projections);
hatI_curr=hatI(:,:,1:N_images);

%%  Global variables

nim=size(hatI_curr,3);
ndef=10; % Number of defocus groups
def1=1;
def2=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter


%% Add CTF and/or envelope
use_CTF=1;
n_im = N_images 
[g_proj_CTF,CTF,index]=  add_CTF_env_v6(hatI_curr(:,:,1:n_im), ndef, def1,def2,B, lambda, use_CTF);

SNR= 1/20;
for count=1:numel(SNR)
    t1_start=tic;
    
    SNR_curr=SNR(count)
    rng(0);
    [noisy_real, noise_v_r]=addnoise_v6(icfft2(g_proj_CTF(:,:,1:n_im)), SNR_curr);
    proj_CTF_noisy=cfft2(noisy_real);
    clear noisy_real
    %% Calculate coeffs using FFBsPCA
    energy_thresh=0.99;
    [ c, R ] = choose_support_v6( proj_CTF_noisy, energy_thresh); %Estimate band limit and compact support size
    c=c*(0.5/floor(size(g_proj_CTF,1)/2)); % Rescaling between 0 and 0.5

    clear  hatI 
    n_r = ceil(4*c*R);
    tic_basis=tic;
    [ basis, sample_points ] = precomp_fb( n_r, R, c );
    timing.basis=toc(tic_basis)
    num_pool=10;
    
    L0=size(g_proj_CTF,1);

    low=0; high=0;
    for nb=1:nbatch
	low=(nb-1)*(K/nbatch)+1;
	high=min(K,low+(K/nbatch)-1);
	sprintf('Saving noisy, demeaned images from %d to %d', low, high)
	curr_batch=proj_CTF_noisy(:,:,low:high);	
	filename=fullfile(fpath_user, sprintf('set%d',nb));
	save(filename,'curr_batch','-v7.3');
    end

    regu=1;
    mean_image_f=mean_LS_batch(CTF, index, regu, K, nbatch, L0, fpath_user);
    err_mean=norm(mean_image_f-mean(hatI_curr,3),'fro')/norm(mean(hatI_curr,3),'fro')

    demean_y_batch(CTF, mean_image_f, index, nbatch, K, fpath_user);
    [ coeff_mean ] = coeff_demean( icfft2(mean_image_f) , R, basis, sample_points, num_pool);
    low=0; high=0;
    for nb=1:nbatch
	low=(nb-1)*(K/nbatch)+1;
	high=min(K,low+(K/nbatch)-1);
	sprintf('Computing coeffs of images from %d to %d', low, high)
	filename=fullfile(fpath_user, sprintf('demean_set%d',nb));
	load(filename);
	[ coeff_pos_k ] = coeff_demean( icfft2(curr_batch) , R, basis, sample_points, num_pool);
	filename=fullfile(fpath_user, sprintf('coeff_set%d',nb));
	save(filename,'coeff_pos_k','-v7.3');
    end

    for k=unique(basis.ang_freqs)'
    sz_k= numel(basis.ang_freqs(basis.ang_freqs==k)); 	
    coeff_pos_k_new=zeros(sz_k,K);
    	for nb=1:nbatch
		low=(nb-1)*(K/nbatch)+1;
		high=min(K,low+(K/nbatch)-1);
		filename=fullfile(fpath_user, sprintf('coeff_set%d',nb));
	        load(filename);
		coeff_pos_k_new(1:sz_k,low:high)=coeff_pos_k{k+1};		
	end
    sprintf('Save coeffs for k=%d',k)
    filename=fullfile(fpath_user, sprintf('coeff_pos%d',k));
    save(filename,'coeff_pos_k_new','-v7.3');
    end
  
    %% CTF in new basis: numerical integration
    [ctf_rad_all]=  calc_CTF_rad(use_CTF, L0, index, ndef, def1,def2,B, lambda, sample_points.r*((floor(L0/2)+1)/0.5));
    jobscript_CCWF_cgshrink_test_batch(index, ctf_rad_all, basis, sample_points,  coeff_mean, noise_v_r, nbatch, nim, fpath_user);

%    energy_coeff=cumsum(cellfun(@norm,denoised_coeff_ccwf));
%    cutoff_coeff=find(energy_coeff/max(energy_coeff)>0.99,1,'first');
%    compressed_den_coeff=denoised_coeff_ccwf(1:cutoff_coeff);
%    sprintf('Total number of coeffs per image after compression is %d' ,sum(cellfun(@numel,compressed_den_coeff))/N_images)

    for nb=1:nbatch
        low=(nb-1)*(K/nbatch)+1;
	high=min(K,low+(K/nbatch)-1);
	for k=unique(basis.ang_freqs)'
		filename=fullfile(fpath_user, sprintf('den_coeff_pos%d',k));
	        load(filename);
		den_coeff_new{k+1}=denoised_coeff(:,low:high);		
	end
    sprintf('Save reshaped batch coeffs for batch%d',nb)
    filename=fullfile(fpath_user, sprintf('batch_den_coeff_pos%d',nb));
    save(filename,'den_coeff_new','-v7.3');
    clear den_coeff_new
    end
  
    [recon] = recon_images_FB_batch(c, R, L0, 1); % Specify batch to reconstruct
%    [mse_ccwf] = calc_MSE_v6(recon,  g_projections(:,:,1:n_im),R)
%    t1_stop=toc(t1_start);
%    mse_ccwf_snr_trend(count)=mse_ccwf;

end 

%% Save results
