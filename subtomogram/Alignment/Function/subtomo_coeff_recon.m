% Input an analytical function in real space, NUFFT it to the Fourier 
% space, sample it at spheres with
% "Lebedev" points, and "Gaussian quadrature" on radial direction. Then
% express the radial component in Fourier Bessel basis, and spherical
% component with spherical harmonic basis.

% %% Gaussian volume to be interpreted
% clear all; close all;
% format long;
% % Gaussian function
% mu = [1 1 0];
% Sigma = [7 8 8];
% % Sigma = [8 1 0; ...
% %          1 7 0; ...
% %          0 0 8]; % to get e-16 in Fourier mvnpdf(pi*[.5 .5 sqrt(2)/2], mu, 1./([8 7 8]))
% 
% % compact support radius to e-15
% thresh = 4; % 4.5 or 5 to attain e-16
% R = round(max(max(Sigma))*thresh);
% 
% % random function not bandlimited
% %F3d = randn(nx,nx,nx)+1i*randn(nx,nx,nx); 
% 
% % sinusoid functions not compactly supported in Fourier
% %F3d = cos(2*X);% + sin(2*Y) + sin(3*Z);
% %% sample volume in real domain and NUFFT transform to Fourier domain
% %nx = round(M^(1/3)/2)*2-1;
% %x1 = linspace(-R,R,nx);
% nx = 2*R+1;
% x1 = -R:1:R;
% [X,Y,Z] = meshgrid(x1,x1,x1);
% 
% F = mvnpdf([X(:) Y(:) Z(:)],mu,Sigma);
% F3d = reshape(F,nx,nx,nx);

function [Error, Time, L] = subtomo_coeff_recon(F3d, R, c, L)
%% Prepare SB basis and Gaussian quadrature points
%c = max(1./Sigma)*thresh; %bandlimit (0, 0.5] has to be < pi
%R = 32;  %Compact support radius/ dimension of volume
nx = 2*R+1;

%c = pi;
%L = 10;
n_r = ceil(2*c*R); %2cR according to Jane, 4cR Nyquist? or sqrt(3)*c according to SB
tic
if ~exist('L','var')
    [ basis, sample_points, L ] = precomp_fb_sphere( n_r, R, c, 2); % 2 checked to third digit
else
    [ basis, sample_points, ~] = precomp_fb_sphere(n_r,R,c,2,L);
end
t_basisprep = toc;
%L = length(basis.Phi_ns); % choosed to match max_angular frequency in the following section
M = (2*L-1)*L*n_r;

%% sample on the sphere in Fourier domain
[thetas,phis] = ssht_sampling(L,'Grid',true);
[x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));
omega = [x y z];

omega_all = omega*sample_points.r(1);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i));
end;

%% forward NUFFT
tic
fg = nufft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,1e-15,nx,nx,nx,F3d);
t_nufftfor = toc;

fg = fg/(2*pi)^3;

display('max fg in outermost sphere',num2str(max(abs(fg(end-L*(2*L-1)+1:end)))));
% fg2 = dirft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,nx,nx,nx,F3d);
% norm(fg-fg2) 1.8357e-12

%% inverse Fourier transform to real domain
%[THETA, FI] = nsht_sampling_points(L);
W = [];% weight for inverse nufft
dphi = 2*pi/(2*L-1);
dtheta = dphi;
for k = 0:n_r-1
    rk = sample_points.r(k+1);
    wk = sample_points.w(k+1);
    w = zeros(L,1);
    for t = 1:L
        w(t) = rk^2*wk*sin(thetas(t,1))*dtheta*dphi; % sin(pi) = 0;
    end
    w(1) = w(1)/2;
    w = repmat(w,2*L-1,1);
    W = cat(1,W,w);
end
%E = sum(W) - 4/3*pi*c^3;

% figure();
% plot(1:L*(2*L-1),W(1:L*(2*L-1)));
% figure();
% plot(1:L,thetas(:,1),1:2*L-1,phis(1,:));

tic
Fg = nufft3d1(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),fg.*W,1,1e-15,nx,nx,nx);
t_nufftback = toc;
% Fg2 = dirft3d1(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),fg.*W,1,nx,nx,nx);
% norm(Fg(:)-Fg2(:)) 8.0168e-14

fac2 = M;%/(2*pi)^3; % scaling to compensate 1/M in nufft1 and volume of input cube
Fg = fac2*reshape(Fg,nx,nx,nx);
err_nufft = norm(complex(F3d(:))-Fg(:))/norm(F3d(:)); %e-3
err_nufft_abs = max(abs(complex(F3d(:))-Fg(:))); %e-6 tried interconversion with filtered, but still same error

N = nx^3;
figure()
subplot(1,3,1);
plot(1:N,real(F3d(:)),1:N,imag(F3d(:)));
title('original volume');
subplot(1,3,2);
plot(1:M,real(fg(:)),1:M,imag(fg(:)));
title('Greengard FT');
subplot(1,3,3);
plot(1:N,real(Fg(:)),1:N,imag(Fg(:)));
title('recon');

%% SH transform on angular direction
flm = zeros(n_r,L^2);
tic
for i = 1:n_r
    temp = (i-1)*L*(2*L-1);
    flm(i,:) = ssht_forward(reshape(fg(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
end
t_shfor = toc;

display('max flm in highest frequency',num2str(max(max(abs(flm(:,(L-1)^2+1:end))))));
%% inverse spherical harmonics transform
tic
F_inv = zeros(M,1);
for i = 1:n_r
    temp = (i-1)*L*(2*L-1);
    F_inv(temp+1:temp+L*(2*L-1)) = reshape(ssht_inverse(flm(i,:),L),M/n_r,1);
end
err_sh = norm(fg-F_inv)/norm(fg); %e-15
err_sh_abs = max(max(abs(fg-F_inv))); %e-17

t_shback = toc;
%% FB transform on radial direction
tic
[ coeff_pos_k ]= FBcoeff_sphere(flm, basis, sample_points);

t_fbfor = toc;
%% inverse FB transform
tic
[ data ]= FBcoeff_sphere_inv(coeff_pos_k, basis, sample_points);
err_fb = norm(data-flm)/norm(flm);   % e-15
[err_fb_abs,~] = max(max(abs(data-flm))); %e-16 be careful if flm is compactly supported!

t_fbback = toc;

% figure();
% plot(1:n_r,real(flm(:,1)),1:n_r,real(data(:,1)),'o');
% figure();
% plot(1:n_r,abs(data(:,1)-flm(:,1)))
% figure();
% plot(1:n_r,sample_points.r);

%% analytical inverse spherical Fourier transform
[ ~, err, t_ift ] = IFT_SB(R, c, L, coeff_pos_k, F3d(:));

Error = [err_sh err_fb err_nufft err(1); err_sh_abs err_fb_abs err_nufft_abs err(2)];
Time = [t_shfor t_fbfor t_nufftfor t_basisprep; t_shback t_fbback t_nufftback t_ift];
end

% %% Error backward from what we can expand
% % FB forward from what we had
% format long;
% tic
% [ coeff_pos_k2 ]= FBcoeff_sphere(data, basis, sample_points);
% toc
% 
% err_fb_31 = 0;
% ind = 0;
% for l = 1:L
%     for m = 1:2*l-1
%         ind = ind+1;
%         err_fb_31 = max(err_fb_31,max(coeff_pos_k{l,ind}-coeff_pos_k2{l,ind}));
%     end
% end
% % e-17
% 
% % FB backward
% [ data2 ]= FBcoeff_sphere_inv(coeff_pos_k2, basis, sample_points);
% err_fb2 = norm(data-data2)/norm(data2);   % e-15
% 
% % SH backward
% F_inv2 = zeros(M,1);
% for i = 1:n_r
%     temp = (i-1)*L*(2*L-1);
%     F_inv2(temp+1:temp+L*(2*L-1)) = reshape(ssht_inverse(data(i,:),L),M/n_r,1);
% end
% % SH forward
% flm2 = zeros(n_r,L^2);
% for i = 1:n_r
%     temp = (i-1)*L*(2*L-1);
%     flm2(i,:) = ssht_forward(reshape(F_inv2(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
% end
% err_sh_abs2 = max(max(abs(data-flm2))); % e-18
% 
%% experiment with pure basis
% tic
% [ coeff_pure ]= FBcoeff_sphere(data, basis, sample_points);
% toc
% 
% ind = 0;
% for l = 1:L
%     for m = 1:2*l-1
%         ind = ind+1;
%         n = length(coeff_pure{l,ind});
%         mag = 0;%max(coeff_pos_k3{l,ind});
%         %mag = 10;
%         coeff_pure{l,ind} = rand(n,1)*mag;
%     end
% end

% % FB backward
% [ data3 ]= FBcoeff_sphere_inv(coeff_pos_k3, basis, sample_points);
% [ coeff_pos_k3_2 ]= FBcoeff_sphere(data3, basis, sample_points);
% err_fb_31 = 0;
% ind = 0;
% for l = 1:L
%     for m = 1:2*l-1
%         ind = ind+1;
%         err_fb_31 = max(err_fb_31,max(coeff_pos_k3{l,ind}-coeff_pos_k3_2{l,ind}));
%     end
% end
% % e-14 mag = max(coeff_pos_k3{l,ind});
% % e-12 mag = 1;
% % e-11 mag = 10;
% % figure()
% % plot(coeff_pos_k3{1,1});
% 
% % SH backward
% F_inv3 = zeros(M,1);
% for i = 1:n_r
%     temp = (i-1)*L*(2*L-1);
%     F_inv3(temp+1:temp+L*(2*L-1)) = reshape(ssht_inverse(data3(i,:),L),M/n_r,1);
% end
% % SH forward
% flm3 = zeros(n_r,L^2);
% for i = 1:n_r
%     temp = (i-1)*L*(2*L-1);
%     flm3(i,:) = ssht_forward(reshape(F_inv3(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
% end
% err_sh_abs3 = max(max(abs(data3-flm3))); % e-17
% 
% % SB forward
% [ coeff_pos_k3_3 ]= FBcoeff_sphere(flm3, basis, sample_points);
% err_fb_32 = 0;
% ind = 0;
% for l = 1:L
%     for m = 1:2*l-1
%         ind = ind+1;
%         err_fb_32 = max(err_fb_32,max(coeff_pos_k3{l,ind}-coeff_pos_k3_3{l,ind}));
%     end
% end
% % e-11