% Generate projections from IP3 molecule, add colored noise and CTF and
% denoise them with CTF after prewhitening.
vol=load('cleanrib.mat');

K=1000;
[projections] = cryo_gen_projections(35,K,1,0);

N_images=K;
g_projections=projections;
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

open_log('test_log.txt')

%% Add CTF and/or envelope
use_CTF=1;
[g_proj_CTF,CTF,index]=  add_CTF_env_v6(hatI_curr, ndef, def1,def2,B, lambda, use_CTF);
clear hatI

%% Loop for different SNR
SNR= [1/10]
for count=1:numel(SNR)
    t1_start=tic;
    
    SNR_curr=SNR(count)
    rng(0);   
    [~, noise, ~, ~] = cryo_addnoise( real(icfft2(g_proj_CTF)),...
        SNR_curr, 'color');
    noisy_real=real(icfft2(g_proj_CTF)) + noise;
    PFDprojs=noisy_real;
    clear noisy_real
    % Normalize images
    %log_message('Normalize background');
    n=size(PFDprojs,1);
    % No need to normalize for synthetic data, needed for real
    % data since exposure may vary across micrographs
    % PFDprojs=cryo_normalize_background(PFDprojs,round(n/2)-10);
    psd = cryo_noise_estimation(PFDprojs);
    %plot(psd(n,:));
    %title('Noise spectrum before prewhitening');
    %Prewhiten
    
    [prewhitened_projs, whiten_filter, nzidx] = cryo_prewhiten(PFDprojs, psd);
    clear PFDprojs
    %psd_white=cryo_noise_estimation(prewhitened_projs);
    %plot(psd_white(n,:));
    %title('Noise spectrum of prewhitened-projections');
    
    proj_CTF_noisy=cfft2(prewhitened_projs);
    [ noise_v_r ] = estimate_noise_real(prewhitened_projs);
    
    % Image is divided by the filter elementwise, so -1
    w_f=cryo_downsample(whiten_filter,size(proj_CTF_noisy,1)).^(-1);
    
    %% Calculate coeffs using FFBsPCA
    n_im = K; % number of denoised images
    energy_thresh=0.99;
    [ c, R ] = choose_support_v6( proj_CTF_noisy, energy_thresh); %Estimate band limit and compact support size
    c=c*(0.5/floor(size(g_proj_CTF,1)/2)); % Rescaling between 0 and 0.5
    
    clear  hatI
    n_r = ceil(4*c*R);
    tic_basis=tic;
    [ basis, sample_points ] = precomp_fb( n_r, R, c );
    timing.basis=toc(tic_basis)
    
    L0=size(g_proj_CTF,1);
    %% "Effective CTF" after the whitening filter
    w_CTF=CTF.*repmat(w_f,1,1,ndef);
    regu=1;
    mean_image_f=mean_LS(CTF, index, proj_CTF_noisy,regu);
    mean_image_f=mean_image_f./w_f;
    err_mean=norm(mean_image_f-mean(hatI_curr,3),'fro')/norm(mean(hatI_curr,3),'fro')
    
    [y_mu] = demean_y_v6(proj_CTF_noisy, w_CTF, mean_image_f, index);
    tic_coeffymu=tic;
    [ coeff_ymu ] = coeff_demean( icfft2(y_mu) , R, basis, sample_points, num_pool );
    timing.coeffymu=toc(tic_coeffymu)
    [ coeff_mean ] = coeff_demean( icfft2(mean_image_f) , R, basis, sample_points, num_pool );
    
    
    %% CTF in new basis: numerical integration
    
    if mod(L0,2)==1
        w_f_rad=interp1([0:floor(L0/2)],w_f(floor(L0/2)+1,floor(L0/2)+1:end),sample_points.r*((floor(L0/2))/0.5) ,'linear');
    else
        w_f_rad=interp1([0:floor(L0/2)],w_f(floor(L0/2),floor(L0/2):end),sample_points.r*((floor(L0/2))/0.5) ,'linear');
    end
    
    
    [ctf_rad_all]=  calc_CTF_rad(use_CTF, L0, index, ndef, def1,def2,B, lambda,sample_points.r*((floor(L0/2))/0.5));
    w_ctf_rad_all=repmat(w_f_rad,1,ndef).*ctf_rad_all;
    
    [ctf_rad_all]=  calc_CTF_rad(use_CTF, L0, index, ndef, def1,def2,B, lambda, sample_points.r*((floor(L0/2)+1)/0.5));
    [ denoised_coeff_ccwf] = jobscript_CCWF_cgshrink_jsb(index, w_f_rad, ctf_rad_all, basis, sample_points,  coeff_mean, coeff_ymu, noise_v_r);
    [recon] = recon_images_FB(c, R, L0, denoised_coeff_ccwf, 1, n_im);
    [mse_ccwf] = calc_MSE_v6(recon,  g_projections(:,:,1:n_im),R)
    
    %% Save results
   
end

