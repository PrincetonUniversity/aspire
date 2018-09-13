function [ images_nufft ] = computeNfftPotts( images, usFftPts, L, points_inside_the_circle)
% This function computes the NFFT of the images at the designated nodes
nImages = size(images,2);
N=2*L;
currImage = zeros(N);
M=size(usFftPts,1);
x = usFftPts.'/pi/2;
plan = nfft_init_guru(2,N,N,M,2*N,2*N,6,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_ESTIMATE);
nfft_set_x(plan,x);
nfft_precompute_psi(plan);

images_nufft = zeros(M, nImages);
for i = 1:nImages    
    currImage(points_inside_the_circle) = images(:,i);
    nfft_set_f_hat(plan,double(currImage(:)));     
    nfft_trafo(plan);
    res = nfft_get_f(plan);
    images_nufft(:,i) = res;  
end    
nfft_finalize(plan);

end

