% Denoise projections of 80s ribosome (with white noise and CTF added).
vol=load('cleanrib.mat');

K=1000;
[projections] = cryo_gen_projections(35,K,1,0);

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
num_pool=num_cores();

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
    
    L0=size(g_proj_CTF,1);

    regu=1;
    mean_image_f=mean_LS(CTF, index, proj_CTF_noisy,regu);
    err_mean=norm(mean_image_f-mean(hatI_curr,3),'fro')/norm(mean(hatI_curr,3),'fro')

    [y_mu] = demean_y_v6(proj_CTF_noisy, CTF, mean_image_f, index);
    [ coeff_ymu ] = coeff_demean( icfft2(y_mu) , R, basis, sample_points, num_pool);
    [ coeff_mean ] = coeff_demean( icfft2(mean_image_f) , R, basis, sample_points, num_pool );
   
    %% CTF in new basis: numerical integration

    [ctf_rad_all]=  calc_CTF_rad(use_CTF, L0, index, ndef, def1,def2,B, lambda, sample_points.r*((floor(L0/2)+1)/0.5));
    [ denoised_coeff_ccwf, skip_flags] = jobscript_CCWF_cgshrink_test(index, ctf_rad_all, basis, sample_points,  coeff_mean, coeff_ymu, noise_v_r);

    energy_coeff=cumsum(cellfun(@norm,denoised_coeff_ccwf));
    cutoff_coeff=find(energy_coeff/max(energy_coeff)>0.99,1,'first');
    compressed_den_coeff=denoised_coeff_ccwf(1:cutoff_coeff);
    sprintf('Total number of coeffs per image after compression is %d' ,sum(cellfun(@numel,compressed_den_coeff))/N_images)


    [recon] = recon_images_FB(c, R, L0, denoised_coeff_ccwf, 1, 1); % Specify range of images to reconstruct
    [mse_ccwf] = calc_MSE_v6(recon,  g_projections(:,:,1:n_im),R)
    t1_stop=toc(t1_start);
    mse_ccwf_snr_trend(count)=mse_ccwf;

end 

%% Save results
