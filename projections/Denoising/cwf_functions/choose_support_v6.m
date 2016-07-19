function [c_limit, R_limit] = choose_support_v6( proj_CTF_noisy, energy_threshold)
%% Determine sizes of the compact support in both real and Fourier space.
% OUTPUTS:
% c_limit: Size of support in Fourier space
% R_limit: Size of support in real space
% We scale the images in real space by L, so that the noise variance in
% both real and Fourier domains is the same.
% Tejal Oct 2015

L=size(proj_CTF_noisy,1);
N=floor(L/2);
P=size(proj_CTF_noisy,3);
[ x, y ] = meshgrid(-N:N, -N:N);
r = sqrt(x.^2 + y.^2);
r_max=N;


img_f=proj_CTF_noisy;
img=(icfft2(img_f))*L;  %% Note: don't use real here because noise variance estimate will then be wrong in line 21, only
%from real part.
mean_data = mean(img, 3);
%remove mean from the data
img = bsxfun(@minus, img, mean_data);

img_corner=reshape(img, L^2, P);
img_corner=img_corner(r>r_max, :);
var_img=var(img_corner(:));

imgf_corner=reshape(img_f, L^2, P);
imgf_corner=imgf_corner(r>r_max, :);
var_imgf=var(imgf_corner(:));

noise_var=min(var_img,var_imgf); %% Note, theoretical img_f and img should give the same variance but there is a small difference, choose the smaller one
% so that you don't get a negative variance or power spectrum in 46,47

variance_map = var(img, [], 3);
%mean 2D variance radial function
radial_var = zeros(N, 1);
for i = 1:N
    radial_var(i) = mean(variance_map(r>=i-1 & r<i));
end;

%for i = 1:size(img,3)
%    img_ps(:, :, i) = abs((img_f(:, :, i))).^2;
%end;
img_ps=abs(img_f).^2;
pspec = mean(img_ps, 3);
radial_pspec = zeros(N, 1);
%compute the radial power spectrum;
for i = 1:N
    radial_pspec(i) = mean(pspec(r>=i-1 & r<i));
end;

%subtract the noise variance
%figure; plot(radial_pspec);
radial_pspec = radial_pspec-noise_var;
radial_var = radial_var-noise_var;
% compute the cumulative variance and power spectrum.
% c = linspace(0, 0.5, N);
% R = 0:N-1;


c = linspace(0,0.5,N)';
R = ((0:N-1))';
cum_pspec = zeros(N, 1);
cum_var = zeros(N, 1);
for i = 1:N
    cum_pspec(i) = sum(radial_pspec(1:i).*c(1:i));
    cum_var(i) = sum(radial_var(1:i).*R(1:i));
end;
cum_pspec = cum_pspec/cum_pspec(end);
cum_var = cum_var/cum_var(end);

%plot(cum_pspec)
%figure; plot(cum_var)
c_limit = c(find(cum_pspec>energy_threshold, 1)-1)*L;
R_limit = R(find(cum_var>energy_threshold, 1) - 1);
