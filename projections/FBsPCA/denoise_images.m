function [ I, Eigim, Freqs, Rad_Freqs, W ] = denoise_images(data, r_max,  K)
%Denoising images: main function
%   Input:  data of size LxLxP
%           r_max, radius of the mask
%           K, number of denoised images wanted 
%   Output: I, denoised images
%           Eigim: eigenimages
%           Freqs: the angular frequencies associated with the eigenimages
%           Rad_Freqs: the associated radial frequencies
%           W: weights for the coefficients (optimal filter)
%Zhizhen Zhao 2012 Sept 17
%set up parameters

P=size(data, 3);
if nargin<3
    K=P;
end;
L=size(data, 1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
%finish set up paramters

%Estimate noise variance from the corners
test=reshape(data, L^2, P);
test=test(r>r_max, :);
noise_variance=var(test(:));
%FBsPCA
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(data, r_max); %Generate steerable PCA basis
[ UU, Freqs, Rad_Freqs, W ] = FBSPCA_MP_rankEst( P, U, D, freqs, rad_freqs, max(noise_variance, D(300))); %Estimate the number of components
[ Coeff ] = WF_FBSPCA( data, Mean, r_max, [0, 0], UU, Freqs, W, 1); %compute the expansion coefficients

%reconstruct the images.
tmp = 2*real(UU(:, Freqs~=0)*Coeff(Freqs~=0, :));
tmp = tmp + UU(:, Freqs==0)*real(Coeff(Freqs==0, :));

I = zeros(L, L, K);
tmp2 = zeros(L);
for i=1:K
    tmp2(r<=r_max)=tmp(:, i);
    I(:, :, i)=tmp2+Mean;
end;

%Eigen images
[~, id]=sort(W, 'descend');
Freqs=Freqs(id);
UU=UU(:, id);
Eigim=zeros(L, L, length(Freqs));
count=1;
for i=1:length(Freqs)
    tmp2(r<=r_max)=real(UU(:, i));
    Eigim(:, :, count)=tmp2;
    count=count+1;
    tmp2(r<=r_max)=imag(UU(:, i));
    Eigim(:, :, count)=tmp2;
end;