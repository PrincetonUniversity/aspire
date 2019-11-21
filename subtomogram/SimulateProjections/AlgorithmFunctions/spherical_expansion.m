function [coeff, basis, sample_points, L] = spherical_expansion(vol, n_r, ss)

% Expand a discrete three-dimensional object on uniform Cartesian grid to 
% spherical Bessel coefficients
%
% INPUT: 
%   vol: volume of size L * L * L, of type double
%   n_r: number of sample points along the radius
%   ss: factor for precision requirement, 3 gives machine precision for 
%       Guassian volumes.
% OUTPUT: 
%   coeff: coefficients with L_max cells, each of size N(l) * (2*l+1), where
%       l is the specific angular frequency and N(l) is the number of
%       corresponding radial basis.
%   basis: spherical Bessel basis with L_max cells, each of size N(l)* 4cR,
%       as they are sampled on Gaussian quadrature points
%   sample_points: 2cL Gaussian quadrature samples and weights evaluated 
%       on [0,c], 

%% precomputation, independent of N
% Prepare SB basis and Gaussian quadrature points
R = floor(size(vol,1)/2);
c = 1/2;

% precompute
[ basis, sample_points, L ] = precomp_fb_sphere(n_r, R, c, ss); 

%% sample on the sphere in Fourier domain
[thetas,phis] = ssht_sampling(L,'Grid',true);
[x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));% thetas(:)-pi/2
omega = [x y z];
omega_all = omega*sample_points.r(1)*2*pi;
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
for i = 1:n_r
    startidx = (i-1)*L*(2*L-1);
    flm1(i,:) = ssht_forward(reshape(fg1(startidx+1:startidx+L*(2*L-1)),L,2*L-1),L);
        %,'Reality',true 
end
clear fg1;
%% FB transform on radial direction
[ coeff ]= FBcoeff_sphere(flm1, basis, sample_points);

end

