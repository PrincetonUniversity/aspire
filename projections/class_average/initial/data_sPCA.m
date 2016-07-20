function [sPCA_data, sPCA_coeff, basis, recon_spca] =  data_sPCA(images, noise_v_r)
% Tejal April 2016

n = size(images, 3);
energy_thresh=0.99;
[ c, R ] = choose_support_v6( cfft2(images), energy_thresh); %Estimate band limit and compact support size
c=c*(0.5/floor(size(images,1)/2)); % Rescaling between 0 and 0.5
c=0.5    
sprintf('hard coded c=%f for testing',c)
n_r = ceil(4*c*R);
tic_basis=tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis=toc(tic_basis)
num_pool=5;

[ timing, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(images, R, noise_v_r, basis, sample_points, num_pool);

sPCA_data.U = U;
sPCA_data.Coeff = cell2mat(sPCA_coeff);
sPCA_data.Mean = mean_coeff;

for i=1:length(sPCA_coeff)
	size_vec(i)=size(sPCA_coeff{i},1);
end


uniq_freq=unique(basis.ang_freqs);
so_far=0;
for i=1:max(uniq_freq)
 if size_vec(i)~=0
   for j=1+so_far:size_vec(i)+so_far
       freqs1(j)=uniq_freq(i);
   end
   so_far=so_far+size_vec(i);
 end
end

sPCA_data.Freqs=freqs1';
sPCA_data.c=c;
sPCA_data.R=R;

L0=size(images,1);
n_max=size(images,3); % Number of images to denoise
%Computes eigen images, need output from IFT_FB.m.
[ fn ] = IFT_FB(R, c);
[~, recon_spca] = denoise_images_analytical(U, fn, mean_coeff, sPCA_coeff, L0, R, n_max);
