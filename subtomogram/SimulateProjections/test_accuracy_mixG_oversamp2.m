addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire');
initpath;
addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master',...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht/'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/nufft3d-1.3.2');
addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'));
addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/nufft2d-1.3.2'));
addpath(genpath('.'));

addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/imrotate3_fast');
addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/spider/bin');

%% test accuracy
m1 = [0,0,-1.04]';
m2 = [0.36, -0.41, 2.56]';
m3 = [-0.54, 0.99, 1.65]';
s1 = [1.9, 1.1, 1.3]';
s2 = [1.4, 0.8, 0.9]';
s3 = [0.3, 0.9, 1.65]';
A = 10;

g = @(x, y, z, sig) A*exp(-((x)/sig(1)).^2/2).* ...
            exp(-((y)/sig(2)).^2/2).* ...
            exp(-((z)/sig(3)).^2/2); % in real space
R0 = 32;

% generate each Gaussian with rotation, translation, and rotation
a_all = [0,0,0];
r_all = eul2rotm(a_all, 'ZYZ');
Cart_orig = @(Cart, m, angle, r) transpose(eul2rotm(angle, 'ZYZ'))*...
    (r'*Cart - repmat(m, 1, size(Cart,2)));

r1 = [0, 65/180*pi, 0];
rot1 = eul2rotm(r1, 'ZYZ');

r2 = [0, 35, -30]./180*pi;
rot2 = eul2rotm(r2, 'ZYZ');

r3 = [162, 54, -37]./180*pi;
rot3 = eul2rotm(r3, 'ZYZ');

bd = 6;
X = linspace(-bd,bd,2*R0+1);
[X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 

% intrinsic with (0, pi/4, 0)
proj_rt = @(x, y) A*(exp((-0.384097 - 0.261152 *x).* x -... 
   0.413223* y.^2).* (1.9426 *erf(13.0713 - 0.12767 *x) + ...
    1.9426 *erf(12.6765 + 0.12767 *x)) +... 
 exp(-0.606377 *x.^2 + ... 
   x.* (2.45316 - 0.124092 *y) + (-0.171666 - ... 
      0.521813 *y).* y) .*(0.108356 *erf( ...
      19.1371 + 0.0628697 *x - 0.357632 *y) + ...
    0.108356 *erf(21.6419 - 0.0628697 *x + 0.357632 *y)) +... 
 exp(-1.20335 *x.^2 + (1.41444 - 1.17474 *y) .*y + ... 
   x .*(0.739245 + 1.16136 *y)) .* (0.430394 *erf( ...
      24.1196 + 1.1984 *x - 1.39898 *y) + ...
    0.430394 *erf(25.5999 - 1.1984 *x + 1.39898 *y)));

fac_com = R0/bd;
G_proj0 = proj_rt(X(:), Y(:)); 
G_proj0 = reshape(G_proj0, 2*R0+1, 2*R0+1);

%% no oversampling
% fac = 1;
% X = linspace(-6,6,2*R/fac+1);
% [X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 
% G_proj2 = proj_rt(X(:), Y(:)); 
% G_proj2 = reshape(G_proj2, 2*R/fac+1, 2*R/fac+1)*fac/fac_com;
% 
% X = linspace(-bd,bd,2*R0/fac+1);
% [X, Y, Z] = ndgrid(X,X,X);% for prolates meshgrid(X, X, X);%ndgrid(X, X);
% Cart = [X(:), Y(:), Z(:)]';
% Cart1 = Cart_orig(Cart, m1, r1, r_all);
% G1 = g(Cart1(1,:), Cart1(2,:), Cart1(3,:), s1);
% Cart2 = Cart_orig(Cart, m2, r2, r_all);
% G2 = g(Cart2(1,:), Cart2(2,:), Cart2(3,:), s2);
% Cart3 = Cart_orig(Cart, m3, r3, r_all);
% G3 = g(Cart3(1,:), Cart3(2,:), Cart3(3,:), s3);
% 
% G = reshape(G1+G2+G3, size(X));    
% 
% R = R0/fac;
% N = 2;
% ns = 1;
% x1 = -R:R;
% nx = length(x1);
% vol = G;
% 
% simulate angles
% NN = N-1; 
% angles = rand(3,NN)*2*pi;
% angles(:,1) = [0, -pi/4, 0];
% angles = [angles zeros(3,1)];%[[0;0;0], [0;3*pi/2;0]];
% 
% compute projections
% [recon, ~, L] = project_volume(vol, angles, R, N, ns, 'real');
% 
% accuracy analysis compare with NUFFT in ASPIRE
% angle2 = [0,pi/4,0];
% rots = eul2rotm(angle2, 'ZYZ');
% projs = cryo_project_np(vol, rots, nx, 'double');
% 
% proj_aspire0 = projs(1:1/fac:end, 1:1/fac:end,1)/fac_com;
% proj_our0 = rot90(flipud(recon(1:1/fac:end, 1:1/fac:end,1)))/fac_com;
% err_our0 = norm(proj_our0 - G_proj0)/norm(G_proj0);%5.654825429459771e-15
% err_aspire0 = norm(proj_aspire0 - G_proj0)/norm(G_proj0);%10^(-8)
% figure()
% subplot(1,2,1)
% imagesc(proj_our0)
% colorbar
% subplot(1,2,2)
% imagesc(G_proj0)
% colorbar
% 
% %  oversample 2x
% fac = 1/3;
% X = linspace(-bd,bd,2*R0/fac+1);
% [X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 
% G_proj1 = proj_rt(X(:), Y(:)); 
% G_proj1 = reshape(G_proj1, 2*R0/fac+1, 2*R0/fac+1);
% 
% X = linspace(-bd,bd,2*R0/fac+1);
% [X, Y, Z] = ndgrid(X,X,X);% for prolates meshgrid(X, X, X);%ndgrid(X, X);
% Cart = [X(:), Y(:), Z(:)]';
% Cart1 = Cart_orig(Cart, m1, r1, r_all);
% G1 = g(Cart1(1,:), Cart1(2,:), Cart1(3,:), s1);
% Cart2 = Cart_orig(Cart, m2, r2, r_all);
% G2 = g(Cart2(1,:), Cart2(2,:), Cart2(3,:), s2);
% Cart3 = Cart_orig(Cart, m3, r3, r_all);
% G3 = g(Cart3(1,:), Cart3(2,:), Cart3(3,:), s3);
% 
% G = reshape(G1+G2+G3, size(X));    
% 
% R = R0/fac;
% N = 2;
% ns = 1;
% x1 = -R:R;
% nx = length(x1);
% vol = G;
% 
% simulate angles
% NN = N-1; 
% angles = rand(3,NN)*2*pi;
% angles(:,1) = [0, -pi/4, 0];
% angles = [angles zeros(3,1)];%[[0;0;0], [0;3*pi/2;0]];
% 
% compute projections
% [recon, ~, L] = project_volume(vol, angles, R, N, ns, 'real');
% 
% accuracy analysis compare with NUFFT in ASPIRE
% angle2 = [0,pi/4,0];
% rots = eul2rotm(angle2, 'ZYZ');
% projs = cryo_project_np(vol, rots, nx, 'double');
% 
% proj_aspire1 = projs(:, :, 1)*fac/fac_com;
% proj_our1 = rot90(flipud(recon(:, :,1)))*fac/fac_com;
% err_our1 = norm(proj_our1 - G_proj1)/norm(G_proj1);%5.654825429459771e-15
% err_aspire1 = norm(proj_aspire1 - G_proj1)/norm(G_proj1);%10^(-8)
% 
% % oversample x3
% fac = 1/5;
% X = linspace(-bd,bd,2*R0/fac+1);
% [X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 
% G_proj2 = proj_rt(X(:), Y(:)); 
% G_proj2 = reshape(G_proj2, 2*R0/fac+1, 2*R0/fac+1);
% 
% X = linspace(-bd,bd,2*R0/fac+1);
% [X, Y, Z] = ndgrid(X,X,X);% for prolates meshgrid(X, X, X);%ndgrid(X, X);
% Cart = [X(:), Y(:), Z(:)]';
% Cart1 = Cart_orig(Cart, m1, r1, r_all);
% G1 = g(Cart1(1,:), Cart1(2,:), Cart1(3,:), s1);
% Cart2 = Cart_orig(Cart, m2, r2, r_all);
% G2 = g(Cart2(1,:), Cart2(2,:), Cart2(3,:), s2);
% Cart3 = Cart_orig(Cart, m3, r3, r_all);
% G3 = g(Cart3(1,:), Cart3(2,:), Cart3(3,:), s3);
% 
% G = reshape(G1+G2+G3, size(X));    
% 
% R = R0/fac;
% N = 2;
% ns = 1;
% x1 = -R:R;
% nx = length(x1);
% 
% simulate angles
% NN = N-1; 
% angles = rand(3,NN)*2*pi;
% angles(:,1) = [0, -pi/4, 0];
% angles = [angles zeros(3,1)];%[[0;0;0], [0;3*pi/2;0]];
% vol = G;
% clear Cart Cart1 Cart2 Cart3 G1 G2 G3 X Y Z G;
% compute projections
% 
% [recon, ~, ~] = project_volume(vol, angles, R, N, ns, 'real');
% 
% accuracy analysis compare with NUFFT in ASPIRE
% angle2 = [0,pi/4,0];
% rots = eul2rotm(angle2, 'ZYZ');
% projs = cryo_project_np(vol, rots, nx, 'double');
% 
% proj_aspire2 = projs(:,:,1)*fac/fac_com;
% proj_our2 = rot90(flipud(recon(:,:,1)))*fac/fac_com;
% err_our2 = norm(proj_our2 - G_proj2)/norm(G_proj2);%5.654825429459771e-15
% err_aspire2 = norm(proj_aspire2 - G_proj2)/norm(G_proj2);%10^(-8)
%% oversample x2
fac = 1/2;
X = linspace(-bd,bd,2*R0/fac+1);
[X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 
G_proj3 = proj_rt(X(:), Y(:)); 
G_proj3 = reshape(G_proj3, 2*R0/fac+1, 2*R0/fac+1);

X = linspace(-bd,bd,2*R0/fac+1);
[X, Y, Z] = ndgrid(X,X,X);% for prolates meshgrid(X, X, X);%ndgrid(X, X);
Cart = [X(:), Y(:), Z(:)]';
Cart1 = Cart_orig(Cart, m1, r1, r_all);
G1 = g(Cart1(1,:), Cart1(2,:), Cart1(3,:), s1);
Cart2 = Cart_orig(Cart, m2, r2, r_all);
G2 = g(Cart2(1,:), Cart2(2,:), Cart2(3,:), s2);
Cart3 = Cart_orig(Cart, m3, r3, r_all);
G3 = g(Cart3(1,:), Cart3(2,:), Cart3(3,:), s3);

G = reshape(G1+G2+G3, size(X));    

R = R0/fac;
N = 2;
ns = 1;
x1 = -R:R;
nx = length(x1);

% simulate angles
NN = N-1; 
angles = rand(3,NN)*2*pi;
angles(:,1) = [0, -pi/4, 0];
angles = [angles zeros(3,1)];%[[0;0;0], [0;3*pi/2;0]];
vol = G;
clear Cart Cart1 Cart2 Cart3 G1 G2 G3 X Y Z G;
% compute projections

[recon, ~, ~] = project_volume(vol, angles, R, N, ns, 'real');

% accuracy analysis compare with NUFFT in ASPIRE
angle2 = [0,pi/4,0];
rots = eul2rotm(angle2, 'ZYZ');
projs = cryo_project_np(vol, rots, nx, 'double');

proj_aspire3 = projs(:,:,1)*fac/fac_com;
proj_our3 = rot90(flipud(recon(:,:,1)))*fac/fac_com;
err_our3 = norm(proj_our3 - G_proj3)/norm(G_proj3);%5.654825429459771e-15
err_aspire3 = norm(proj_aspire3 - G_proj3)/norm(G_proj3);%10^(-8)

%% oversample x4

fac = 1/4;
X = linspace(-bd,bd,2*R0/fac+1);
[X, Y] = ndgrid(X, X);% for prolates meshgrid(X, X, X);%ndgrid(X, X); 
G_proj4 = proj_rt(X(:), Y(:)); 
G_proj4 = reshape(G_proj4, 2*R0/fac+1, 2*R0/fac+1);

X = linspace(-bd,bd,2*R0/fac+1);
[X, Y, Z] = ndgrid(X,X,X);% for prolates meshgrid(X, X, X);%ndgrid(X, X);
Cart = [X(:), Y(:), Z(:)]';
Cart1 = Cart_orig(Cart, m1, r1, r_all);
G1 = g(Cart1(1,:), Cart1(2,:), Cart1(3,:), s1);
Cart2 = Cart_orig(Cart, m2, r2, r_all);
G2 = g(Cart2(1,:), Cart2(2,:), Cart2(3,:), s2);
Cart3 = Cart_orig(Cart, m3, r3, r_all);
G3 = g(Cart3(1,:), Cart3(2,:), Cart3(3,:), s3);

G = reshape(G1+G2+G3, size(X));    

R = R0/fac;
N = 2;
ns = 1;
x1 = -R:R;
nx = length(x1);

% simulate angles
NN = N-1; 
angles = rand(3,NN)*2*pi;
angles(:,1) = [0, -pi/4, 0];
angles = [angles zeros(3,1)];%[[0;0;0], [0;3*pi/2;0]];
vol = G;
clear Cart Cart1 Cart2 Cart3 G1 G2 G3 X Y Z G;
% compute projections

[recon, ~, ~] = project_volume(vol, angles, R, N, ns, 'real');

% accuracy analysis compare with NUFFT in ASPIRE
angle2 = [0,pi/4,0];
rots = eul2rotm(angle2, 'ZYZ');
projs = cryo_project_np(vol, rots, nx, 'double');

proj_aspire4 = projs(:,:,1)*fac/fac_com;
proj_our4 = rot90(flipud(recon(:,:,1)))*fac/fac_com;
err_our4 = norm(proj_our4 - G_proj4)/norm(G_proj4);%5.654825429459771e-15
err_aspire4 = norm(proj_aspire4 - G_proj4)/norm(G_proj4);%10^(-8)


save('mixG6_24.mat');
fprintf('done');
%% print tables for the errors

% err_aspire = zeros(5,3);
% err_our = zeros(5,3);
% load('mixG6_135.mat');
% load('mixG6_24.mat');
% err_aspire(1,1) = err_aspire0;
% err_aspire(3,1) = err_aspire1;
% err_aspire(5,1) = err_aspire2;
% err_aspire(2,1) = err_aspire3;
% err_aspire(4,1) = err_aspire4;
% err_our(1,1) = err_our0;
% err_our(3,1) = err_our1;
% err_our(5,1) = err_our2;
% err_our(2,1) = err_our3;
% err_our(4,1) = err_our4;
% 
% load('mixG16_135.mat');
% load('mixG16_24.mat');
% err_aspire(1,2) = err_aspire0;
% err_aspire(3,2) = err_aspire1;
% err_aspire(5,2) = err_aspire2;
% err_aspire(2,2) = err_aspire3;
% err_aspire(4,2) = err_aspire4;
% err_our(1,2) = err_our0;
% err_our(3,2) = err_our1;
% err_our(5,2) = err_our2;
% err_our(2,2) = err_our3;
% err_our(4,2) = err_our4;
% 
% 
% load('mixG26_12345.mat');
% err_aspire(1,3) = err_aspire0;
% err_aspire(3,3) = err_aspire1;
% err_aspire(5,3) = err_aspire2;
% err_aspire(2,3) = err_aspire3;
% err_aspire(4,3) = err_aspire4;
% err_our(1,3) = err_our0;
% err_our(3,3) = err_our1;
% err_our(5,3) = err_our2;
% err_our(2,3) = err_our3;
% err_our(4,3) = err_our4;
% 
% format shortEng
% format compact
% err_table = array2table(err_aspire, 'RowNames',...
%     {'x1','x2','x3',...
%     'x4','x5'}, ...
%     'VariableNames',{'range6','range16','range26'});
% display(err_table);
% format
% writetable(err_table,'err_mixGaussian_over.csv','Delimiter',',','QuoteStrings',true,...
%     'WriteRowNames',true)
% err_table = array2table(err_our, 'RowNames',...
%     {'x1','x2','x3',...
%     'x4','x5'}, ...
%     'VariableNames',{'range6','range16','range26'});
% display(err_table);
% format
% writetable(err_table,'err_mixGaussian_over_our.csv','Delimiter',',','QuoteStrings',true,...
%     'WriteRowNames',true)

%%
%load('mixG16_24.mat');
% row = 5; col = 3;
% proj = cell(col,row);
% proj{1,1} = G_proj0;
% proj{1,2} = proj_our0;
% proj{1,4} = proj_aspire0;
% proj{2,1} = G_proj1;
% proj{2,2} = proj_our1;
% proj{2,4} = proj_aspire1;
% proj{3,1} = G_proj2;
% proj{3,2} = proj_our2;
% proj{3,4} = proj_aspire2;
% for r = 1:col
%     proj{r,3} = abs(proj{r,2} - proj{r,1})/norm(proj{r,1});
%     proj{r,5} = abs(proj{r,4} - proj{r,1})/norm(proj{r,1});
% end
% 
% 
% names = { 'integration (0,pi/4,0) x1','integration (0,pi/4,0) x2.5',...
%     'integration (0,pi/4,0) x5',...
%     'spherical expansion','spherical expansion','spherical expansion',...
%     'error of spherical expansion','error of spherical expansion',...
%     'error of spherical expansion',...
%     'ASPIRE','ASPIRE','ASPIRE',...
%     'error of ASPIRE','error of ASPIRE','error of ASPIRE' };
% 
% figure()
% count = 0;
% for i = 1:row
%     for j = 1:col
%         count = count+1;
%         subplot(row,col,count)
%         imagesc(proj{j,i});
%         colormap(flipud(gray)); 
%         title(names(count),'FontSize',8)
%         colorbar; axis square;
%         ax = gca;
%         ax.FontSize = 8; 
%         
%     end
% end
% suptitle('-16<x,y,z<16, L0 = 32')
% % print('-bestfit','mixG_over3','-dpdf')
% %print('-fillpage','mixG_over6','-dpdf')