%Denoise example code, input image
var_hat = 1; %noise variance, change accordingly
n_im = 3; % number of denoised images
energy_thresh=0.99;
data=ones(10,10,n_im)+0.1*randn(10,10,n_im);
[ c, R ] = avg_pspec(data, var_hat, energy_thresh); %Estimate band limit and compact support size
num_pool=2;
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_fb( n_r, R, c );
[ timing, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(data, R, var_hat, basis, sample_points, num_pool); % Compute Fast steerable PCA.
[ fn ] = IFT_FB(R, c); % compute inverse fourier transform of the Fourier-Bessel functions
L0=size(data,1);
[ mean_image, denoised ] = denoise_images_analytical(U, fn, mean_coeff, sPCA_coeff, L0, R, n_im); %generate denoised images.


