function [coeff, basis, sample_points, L] = spherical_expansion_ss(vol, n_r, ss, maxL)
% expand a discrete three-dimensional object on Cartesian grid to spherical
% Bessel coefficients
%
% INPUT: 
%   vol: volume of size L * L * L, of type double
% OUTPUT: 
%   coeff: coefficients with L_max cells, each of size N(l) * (2*l+1), where
%       l is the specific angular frequency and N(l) is the number of
%       corresponding radial basis.
%   basis: 
%   sample_points:

%% precomputation, independent of N
% Prepare SB basis and Gaussian quadrature points
R = floor(size(vol,1)/2);
c = 1/2;

% precompute
if ~exist('maxL','var')
    [ basis, sample_points, L ] = precomp_fb_sphere(n_r, R, c, ss); 
else
    [ basis, sample_points, ~ ] = precomp_fb_sphere(n_r, R, c, ss, maxL);
    L = maxL;
end

%% sample on the sphere in Fourier domain
[thetas,phis] = ssht_sampling(L,'Method', 'MWSS', 'Grid',true);
[x, y, z] = sph2cart(phis(:,1:L),thetas(:,1:L)-pi/2,ones((L+1),L));% thetas(:)-pi/2
omega = [x(:) y(:) z(:)];
omega_all = omega*sample_points.r(1)*(-2*pi);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i)*2*pi);
end

% %% sample on the sphere in Fourier domain
% [x, y, z] = nsht_plot_sampling(L);
% omega = [x(:) y(:) z(:)]*2*pi;
% 
% omega_all = omega*sample_points.r(1);
% for i = 2:n_r
%     omega_all = cat(1,omega_all,omega*sample_points.r(i));
% end

fout = omega_all(:,1)>=pi | omega_all(:,2)>=pi | omega_all(:,3)>=pi;
clear x y z;

%% forward NUFFT
set_nufft_libraries('finufft');
fg1 = nufft3(vol, omega_all');
fg1(fout) = 0;
clear omega_all fout;
%E = zeros(2,n_r);
% fg1_new = zeros((L+1)*2*L*n_r,1);
% for i = 1:n_r
% 
% startidx = (i-1)*(L+1)*2*L;
% sphere = reshape(fg1(startidx+1:startidx+(L+1)*L), L+1, L);
% fg1_new(startidx+1:startidx+(L+1)*2*L) = reshape(cat(2,sphere,flipud(conj(sphere))), (L+1)*L*2, 1);
% 
% startidx1 = (i-1)*(L+1)*(L);
% sphere2 = reshape(fg1_all(startidx+1:startidx+(L+1)*2*L), L+1, 2*L);
% E(1,i) = norm(imag(flipud(sphere(:,1:L))) + imag(sphere2(:,L+1:end)));
% E(2,i) = norm(real(flipud(sphere(:,1:L))) - real(sphere2(:,L+1:end)));
% 
% end
% %% SH transform on angular direction
% flm1 = zeros(n_r,L^2);
% for i = 1:n_r
%     startidx = (i-1)*L^2;%(i-1)*L*(2*L-1);
%     flm1(i,:) = nsht_forward(fg1(startidx+1:startidx+L^2), L);
%     %ssht_forward(reshape(fg1(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
%         %,'Reality',true 
% end
% clear fg1;
%% SH transform on angular direction
flm1 = zeros(n_r,L^2);
fg1 = reshape(fg1, (L+1)*L, n_r);
for i = 1:n_r %par
    sphere = reshape(fg1(:,i), L+1, L);
    flm1(i,:) = ssht_forward(cat(2,sphere,flipud(conj(sphere))),L,'Method','MWSS');
        %,'Reality',true 
end
clear fg1;
%% FB transform on radial direction
[ coeff ]= FBcoeff_sphere(flm1, basis, sample_points);

end

