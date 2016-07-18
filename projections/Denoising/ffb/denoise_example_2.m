%Denoise example code, input image
var_hat = 1; %noise variance, change accordingly
n_im = 100; % number of denoised images
data=ones(10,10,n_im)+0.1*randn(10,10,n_im);
energy_thresh=0.99;
[ c, R ] = avg_pspec(data, var_hat, energy_thresh); %Estimate band limit and compact support size
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_fb( n_r, R, c );

