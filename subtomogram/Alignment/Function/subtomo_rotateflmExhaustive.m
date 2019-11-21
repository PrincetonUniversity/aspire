function [E, T] = subtomo_rotateflmExhaustive(alpha, beta, gamma);
% Expand a volume and its rotated volume in Fourier domain
% find the relative rotation by exhaustive search

% cd /u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE %on server

% Gaussian volume to be interpreted
clear all; close all;
format long;
% Gaussian function
mu1 = [0 1 0]; mu1_y = [0 1 0]; mu1_z = [0 0 0];
mu2 = [0 0 1]; mu2_y = [0 0 1]; mu2_z = [0 0 0];
Sigma1 = [7 8 6]; Sigma1_y = [8 7 8]; Sigma1_z = [8 7.3 8];
Sigma2 = [7 6 8]; Sigma2_y = [8 8 7]; Sigma2_z = [8 8 7.3];

% compact support radius to e-15
thresh = 4; % 4.5 or 5 to attain e-16
R = round(max(max(Sigma1))*thresh);

% sample volume in real domain and NUFFT transform to Fourier domain
nx = 2*R+1;
x1 = -R:1:R;
[Y,X,Z] = meshgrid(x1,x1,x1);

F1 = mvnpdf([X(:) Y(:) Z(:)],mu1,Sigma1);
%F1_y = mvnpdf([X(:) Y(:) Z(:)],mu1_y,Sigma1_y);
%F1_z = mvnpdf([X(:) Y(:) Z(:)],mu1_z,Sigma1_z);
F2 = mvnpdf([X(:) Y(:) Z(:)],mu2,Sigma2);
%F2_y = mvnpdf([X(:) Y(:) Z(:)],mu2_y,Sigma2_y);
%F2_z = mvnpdf([X(:) Y(:) Z(:)],mu2_z,Sigma2_z);
F3d1 = reshape(F1,nx,nx,nx); %+F1_z+F1_y
F3d2 = reshape(F2,nx,nx,nx); %+F2_z+F2_y

% ribosome
%load('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/projections/simulation/maps/cleanrib.mat')


%% Prepare SB basis and Gaussian quadrature points
tic
c = pi;
L = 20;
n_r = ceil(2*c*R); %2cR according to Jane, 4cR Nyquist? or sqrt(3)*c according to SB
%[ basis, sample_points, L ] = precomp_fb_sphere( n_r, R, c, 2);
[ basis, sample_points ] = precomp_fb_sphere( n_r, R, c, 2, L); % 2 checked to third digit
M = (2*L-1)*L*n_r;
t_precomp_fb = toc;

% sample on the sphere in Fourier domain
tic
[thetas,phis] = ssht_sampling(L,'Grid',true);
[x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));% thetas(:)-pi/2
omega = [x y z];

omega_all = omega*sample_points.r(1);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i));
end;
t_precomp_sh = toc;

%% forward NUFFT
fg1 = nufft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,1e-15,nx,nx,nx,F3d1);
tic
fg2 = nufft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,1e-15,nx,nx,nx,F3d2);
fg1 = fg1/(2*c)^3;
fg2 = fg2/(2*c)^3;
t_nufftfor = toc;

display('vol1 max fg in outermost sphere',num2str(max(abs(fg1(end-L*(2*L-1)+1:end)))));
display('vol2 max fg in outermost sphere',num2str(max(abs(fg2(end-L*(2*L-1)+1:end)))));
clear F3d1 F3d2;
%% SH transform on angular direction
flm1 = zeros(n_r,L^2);
for i = 1:n_r
    temp = (i-1)*L*(2*L-1);
    flm1(i,:) = ssht_forward(reshape(fg1(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
end

tic
flm2 = zeros(n_r,L^2);
for i = 1:n_r
    temp = (i-1)*L*(2*L-1);
    flm2(i,:) = ssht_forward(reshape(fg2(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
end
t_shfor = toc;
display('max flm in highest frequency, vol1',num2str(max(max(abs(flm1(:,(L-1)^2+1:end))))));

display('max flm in highest frequency, vol2',num2str(max(max(abs(flm2(:,(L-1)^2+1:end))))));

%% rotate a single sphere in SH 
% right_hand viewing from positive axis!extrinsic rotation zyz!! rotate the
% volume instead of the axes. and notice there is an additiona pi rotation
% in x on the complex part.
% For y to z: alpha = pi/2;beta = pi/2;gamma = -pi/2;
alpha = pi/2; 
beta = pi/2;
gamma = -pi/2;

% take one of the spheres to visualize
sph_ind = 1;
flm_unrot = flm1(sph_ind,:)';
flm_rotated = flm2(sph_ind,:)';

%define the rotation matrices succesively for plotting
R1 = euler2rotationMatrix(alpha, 0, 0, 'zyz');
R2 = euler2rotationMatrix(alpha, beta, 0, 'zyz');
R3 = euler2rotationMatrix(alpha, beta, gamma, 'zyz');

% rotation matrices using the recursive method for complex SHs
R_SH1 = getSHrotMtx(R1, L-1, 'complex');
R_SH2 = getSHrotMtx(R2, L-1, 'complex');
R_SH3 = getSHrotMtx(R3, L-1, 'complex');

% rotated coefficients of complex function
flm_rot1 = R_SH1*flm_unrot;
flm_rot2 = R_SH2*flm_unrot;
flm_rot3 = R_SH3*flm_unrot;

clear R_SH1 R_SH2;
 
% % check with vol1 to see if the angles are right
max(flm_rot3-flm_rotated) 

%Compute rotated function on the sphere.
f_rot1 = ssht_inverse(complex(real(flm_rot1), imag(flm_rot1)), L);
f_rot2 = ssht_inverse(complex(real(flm_rot2), imag(flm_rot2)), L);
f_rot3 = ssht_inverse(complex(real(flm_rot3), imag(flm_rot3)), L);


%% rotate the whole volume and find relative orientation with correlation
flm_rot = R_SH3*flm1'; % transpose so that each column represent a sphere
flm_rot = flm_rot';

load('SO3_10.mat');
%load('Eul_20.mat'); % Eul_20 is not accurate, or it's the extrinsic angle
rot = SO3;
num_rot = size(rot,3);
cor_sph = zeros(1,num_rot);
l2_sph = zeros(1,num_rot);

%% for a single sphere
l = 5;
max(flm_unrot(l^2+1:end))

% for volume 1 and 2
tic
for i = 1:num_rot
    l2_sph(i) = norm(flm2(1,1:l^2)' - getSHrotMtx(rot(:,:,i),l-1,'complex')*flm1(1,1:l^2)');
end
t_l2_singlesph = toc;
[m_l2,I_l2_2] = min(l2_sph); % good
rot_l2 = rot(:,:,I_l2_2);
% ang_l2 = Eul_20(I_l2_2,:)/pi*180;

% check if they are equal
% sum_l2 = (norm(flm_rot1(1:l^2) - getSHrotMtx(R1,l-1,'complex')*flm_unrot(1:l^2)))^2;
% sum_cor = norm(flm_rot1(1:l^2))^2+norm(flm_unrot(1:l^2))^2 ...
%  -2*real(conj(flm_rot1(1:l^2)')*(getSHrotMtx(R1,l-1,'complex')*flm_unrot(1:l^2)));
% % good now! there has to be no shift
% display(sum_l2-sum_cor);

% use correlation
% tic
% for i = 1:num_rot
%     cor_sph(i) = flm_unrot(1:l^2)'*getSHrotMtx(rot(:,:,i),l-1,'complex')*conj(flm_rot1(1:l^2)); % worked for imag part, why?
% end
% toc
% [m_cor,I_cor] = max(real(cor_sph));    
% ang_cor = Eul_20(I_cor,:)/pi*180;

tic
for i = 1:num_rot
    cor_sph(i) = real(conj(flm2(1:1:l^2))*(getSHrotMtx(rot(:,:,i),l-1,'complex')*flm1(1,1:l^2)'));
end
t_cor_singlesph = toc;

%[value, index] = sort(cor_sph,'descend');
[~,I_cor] = max(cor_sph);    
%ang_cor1 = Eul_20(I_cor,:)/pi*180;
rot_cor1 = rot(:,:,I_cor);

%% for the whole volume
w_sph = repmat(sample_points.r.^2.*sample_points.w,1,l^2);

tic
for i = 1:num_rot
    cor_sph(i) = trace(conj(flm2(:,1:l^2))*(getSHrotMtx(rot(:,:,i),l-1,'complex'))*(w_sph.*flm1(:,1:l^2))');
end
t_cor_search = toc; % checked but to some other rotation??

[m_cor,I_cor] = max(real(cor_sph));    
rot_cor = rot(:,:,I_cor); % recovered R3'
ang_cor = Eul_20(I_cor,:)/pi*180;

%% FB transform on radial direction
tic
[ coeff1 ]= FBcoeff_sphere(flm1, basis, sample_points);

[ coeff2 ] = FBcoeff_sphere(flm2, basis, sample_points);
t_fbfor = toc;

[coeff_rot] = FBcoeff_sphere(flm_rot, basis, sample_points);

C1 = zeros(length(coeff1{1,1}),l^2);
C2 = C1;
for i = 1:l
    for j = 1:2*i-1
        len = length(coeff1{i,(i-1)^2+j});
        C1(1:len,(i-1)^2+j) = coeff1{i,(i-1)^2+j};
        C2(1:len,(i-1)^2+j) = coeff2{i,(i-1)^2+j};
    end
end

cor_coeff = zeros(1,num_rot);
tic
for i = 1:num_rot
    cor_coeff(i) = real(trace(conj(C2(:,1:l^2))*(getSHrotMtx(rot(:,:,i),l-1,'complex'))*(C1(:,1:l^2))'));
end
t_coeff_search = toc; % checked but to some other rotation??

[m_coeff,I_coeff] = max(real(cor_coeff));    
rot_coeff = rot(:,:,I_coeff);
%ang_coeff = Eul_20(I_coeff,:)/pi*180;


% [coeff_rots] = FBcoeff_sphere(flm_rots, basis, sample_points);

% %% inverse FB transform
% tic
% [ data ]= FBcoeff_sphere_inv(coeff_pos_k, basis, sample_points);
% err_fb = norm(data-flm)/norm(flm);   % e-15
% [err_fb_abs,I] = max(max(abs(data-flm))); %e-16 be careful if flm is compactly supported!
% 
% t_fbback = toc;

%% rotate with the estimated rotation to check its accuracy.
R_SHest = getSHrotMtx(rot_coeff,L-1,'complex');
flm_rotest = R_SHest*flm_unrot;
flm_rotest = flm_rotest';
f_rotest = ssht_inverse(complex(real(flm_rotest),imag(flm_rotest)),L); 

%%
temp = (sph_ind-1)*L*(2*L-1);
f1 = reshape(fg1(temp+1:temp+L*(2*L-1)),L,2*L-1);
f2 = reshape(fg2(temp+1:temp+L*(2*L-1)),L,2*L-1);
type = 'colour';
figure(1);
subplot(2,5,1)
ssht_plot_sphere(real(f1), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
title('1 sph unrot in Fourier');
subplot(2,5,6)
ssht_plot_sphere(imag(f1), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
% Plot rotated functions on sphere.
subplot(2,5,2)
ssht_plot_sphere(real(f_rot1), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
title(['alpha=' num2str(alpha) ' rot']);
subplot(2,5,7)
ssht_plot_sphere(imag(f_rot1), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
subplot(2,5,3)
ssht_plot_sphere(real(f_rot2), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
title(['and beta=' num2str(beta) ' rot']);
subplot(2,5,8)
ssht_plot_sphere(imag(f_rot2), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
subplot(2,5,4)
ssht_plot_sphere(real(f_rot3), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
title(['and gamma=' num2str(gamma) ' rot']);
subplot(2,5,9)
ssht_plot_sphere(imag(f_rot3), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
subplot(2,5,5)
ssht_plot_sphere(real(f_rotest), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
title('1 sph rotated with est');
subplot(2,5,10)
ssht_plot_sphere(imag(f_rotest), L, thetas, phis, 'Type', type, ...
    'ColourBar', false);
hold off;

% rot_coeff is the identical rotation! why alpha seems to be inaccurate?

%% analytical inverse spherical Fourier transform
% [ vol, err, t_ift ] = IFT_SB(R, c, L, coeff1, F1+F1_y);
% [ vol_rot, ~, t_ift_rot] = IFT_SB(R,c,L,coeff_rot);
% [ vol_rots, ~, t_ift_rots] = IFT_SB(R,c,L,coeff_rots);
tic
[ basis2 ] = precomp_IFT_SB(R, c, L);
t_precomp_ift = toc;
tic
[ vol1, err1, t_ift1 ] = coeff_IFT_SB(basis2,coeff1, F1);%+F1_y+F1_z
t_ift = toc;
[ vol_rot, err_rot, t_ift2 ] = coeff_IFT_SB(basis2,coeff_rot, F2);%+F2_y+F2_z

end

%% check sig_rot = R3*Sigma'; works now
% F_r = F2+F2_y;
% [phi, theta, r] = cart2sph(X,Y,Z); % phi is the longitude, theta is the latitude
% r = r(:);
% phi = phi(:);
% theta = theta(:);
% ball = find(r<=R);
% F_r = F_r(ball);
% theta = theta(ball);
% phi = phi(ball);
% %err_rot = max(abs(vol_rots-vol_rot)); %e-17
% err_rot_rel = norm(abs(vol_rot-F_r))/norm(F_r); %e-17
% err_rot_abs = max(abs(vol_rot-F_r)); %e-14
% origin = (length(ball)+1)/2;
% ind = origin:origin+80;
% figure()
% plot(ind,vol_rots(ind),'o-',ind,F_r(ind),'o-');
% plot(theta,vol_rot-F_r);
% plot(phi,vol_rot-F_r);


