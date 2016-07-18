%Denoise example code, input image
var_hat = 1; %noise variance, change accordingly
n_im = 100; % number of denoised images
[ c, R ] = avg_pspec(data, var_hat); %Estimate band limit and compact support size
data=ones(10,10,n_im)+0.1*randn(10,10,n_im);
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_fb( n_r, R, c );

