% Input an analytical function in real space, NUFFT it to the Fourier 
% space, sample it at spheres with
% "Lebedev" points, and "Gaussian quadrature" on radial direction. Then
% express the radial component in Fourier Bessel basis, and spherical
% component with spherical harmonic basis.

function [E,T,L] = subtomo_nufft_sh_FB(R, fac, c)
%% Prepare SB basis and Gaussian quadrature points
E = zeros(2,4);
T = zeros(1,7);
tic
%c = pi; %bandlimit (0, 0.5]
%R = 4;  %Compact support radius/ dimension of volume
n_r = ceil(4*c*R);
[ basis, sample_points ] = precomp_fb_sphere( n_r, R, c, 2); % checked to third digit

%% sample on the sphere in Fourier domain
L = length(basis.Phi_ns); % choosed to match max_angular frequency in the following section
M = L^2*n_r;
[x, y, z] = nsht_plot_sampling(L);
omega = [x(:) y(:) z(:)];

omega_all = omega*sample_points.r(1);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i));
end;
T(1) = toc;

%% sample volume in real domain and NUFFT transform to Fourier domain
nx = round(M^(1/3)/2)*2-1;
x1 = linspace(-pi,pi,nx);
[X,Y,Z] = meshgrid(x1,x1,x1);

% Gaussian function
%fac = 1/4;
mu = [0 0 0];
Sigma = [1 1.2 1.3]*fac;
F = mvnpdf([X(:) Y(:) Z(:)],mu,Sigma);
F3d = reshape(F,nx,nx,nx);

fg = nufft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,1e-15,nx,nx,nx,F3d);
T(2) = toc;

% fg2 = dirft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,nx,nx,nx,F3d);
% norm(fg-fg2) 1.8357e-12

%% inverse Fourier transform to real domain
[THETA, ~] = nsht_sampling_points(L);
W = zeros(M,1);% weight for inverse nufft
dtheta = 2*pi/(2*L-1);
for k = 0:n_r-1
    rk = sample_points.r(k+1);
    wk = sample_points.w(k+1);
    for l = 0:L-1
        dphi = 2*pi/(2*l+1);
        wl = rk^2*wk*sin(THETA(l+1)) *dtheta*dphi;
        W(k*L^2 + l^2 +1: k*L^2 + (l+1)^2) = wl*ones(2*l+1,1);
    end
end
E(1,1) = abs(sum(W) - 4/3*pi*c^3);

Fg = nufft3d1(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),fg.*W,1,1e-15,nx,nx,nx);

% Fg2 = dirft3d1(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),fg.*W,1,nx,nx,nx);
% norm(Fg(:)-Fg2(:)) 8.0168e-14

fac2 = M/(2*pi)^3; % scaling to compensate 1/M in nufft1 and volume of input cube
Fg = fac2*reshape(Fg,nx,nx,nx);
err_nufft = norm(complex(F3d(:))-Fg(:))/norm(F3d(:)); %e-4
err_nufft_abs = norm(complex(F3d(:))-Fg(:));


% N = nx^3;
% figure()
% subplot(1,3,1);
% plot(1:N,real(F(:)),1:N,imag(F(:)));
% title('original volume');
% subplot(1,3,2);
% plot(1:M,real(fg(:)),1:M,imag(fg(:)));
% title('Greengard FT');
% subplot(1,3,3);
% plot(1:N,real(Fg(:)),1:N,imag(Fg(:)));
% title('recon');
E(1,2) = err_nufft;
E(2,2) = err_nufft_abs;
T(3) = toc;
%% SH transform on angular direction
% flm = zeros(n_r,L^2);
% for i = 1:n_r
%     temp = (i-1)*L^2;
%     flm(i,:) = nsht_forward(fg(temp+1:temp+L^2)',L);
% end
% 
% T(4) = toc;
% %% inverse spherical harmonics transform
% F_inv = zeros(L^2*n_r,1);
% for i = 1:n_r
%     temp = (i-1)*L^2;
%     F_inv(temp+1:temp+L^2) = nsht_inverse(flm(i,:),L)';
% end
% err_sh = norm(fg-F_inv)/norm(fg); %e-16
% err_sh_abs = norm(fg-F_inv); %e-12
% 
% E(1,3) = err_sh;
% E(2,3) = err_sh_abs;
% T(5) = toc;
% %% FB transform on radial direction
% [ coeff_pos_k ]= FBcoeff_sphere(flm, basis, sample_points);
% 
% T(6) = toc;
% %% inverse FB transform
% [ data ]= FBcoeff_sphere_inv(coeff_pos_k, basis, sample_points);
% err_fb = norm(data-flm)/norm(flm);   % e-7
% err_fb_abs = norm(data-flm); %e-5
% 
% E(1,4) = err_fb;
% E(2,4) = err_fb_abs;
% T(7) = toc;

end
