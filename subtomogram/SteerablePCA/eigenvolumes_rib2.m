% 0. generate random SB coefficients
% 3. compute its steerable PCA to get eigenvolumes
% 4. compute N rotated version of it
% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
% 6. compare eigenvectors and eigenvalues from classical and steerable PCA

% Dec 1, 2017

addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner');

% ribosome in real domain
load('cleanrib.mat');
R = 21;
F3d = real(volref(33-R:33+R,33-R:33+R,33-R:33+R));
nx = 2*R+1;
x1 = -R:1:R;
[X,Y,Z] = meshgrid(x1,x1,x1);
% Fourier domain
c = pi;
L = R+1;
n_r = ceil(2*c*R); %2cR according to Jane, 4cR Nyquist? or sqrt(3)*c according to SB
M = (2*L-1)*L*n_r;
[thetas,phis] = ssht_sampling(L,'Grid',true);
[x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));% thetas(:)-pi/2
omega = [x y z];
[ basis, sample_points ] = precomp_fb_sphere( n_r, R, c, 2, L); 
omega_all = omega*sample_points.r(1);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i));
end;
fg1 = nufft3d2(M,omega_all(:,1),omega_all(:,2),omega_all(:,3),-1,1e-15,nx,nx,nx,F3d);
fg1 = fg1/(2*c)^3;
display('vol1 max fg in outermost sphere',num2str(max(abs(fg1(end-L*(2*L-1)+1:end)))));

% visualize volume in real space
figure()
subplot(1,3,1)
threshF = quantile(F3d(:),0.1);
Fv = isosurface(X,Y,Z,F3d,threshF);
p = patch(Fv);
%isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('clean E coli ribosome in real space');
% visualize in Fourier space
subplot(1,3,2)
fReal = real(fg1);
fR = scatteredInterpolant(omega_all(:,1),omega_all(:,2),omega_all(:,3),fReal);
[xq,yq,zq] = meshgrid(-pi:0.1:pi);
vq = fR(xq,yq,zq);
threshfR = quantile(fReal,0.1);
fv = isosurface(xq,yq,zq,vq,threshfR);
p = patch(fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('real part in Fourier space');

subplot(1,3,3)
fComp = imag(fg1);
fC = scatteredInterpolant(omega_all(:,1),omega_all(:,2),omega_all(:,3),fComp);
vq = fC(xq,yq,zq);
threshfC = quantile(fComp,0.1);
fv = isosurface(xq,yq,zq,vq,threshfC);
p = patch(fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('complex part in Fourier space');

% convert to SB coefficients
flm1 = zeros(n_r,L^2);
tic
for i = 1:n_r
    temp = (i-1)*L*(2*L-1);
    flm1(i,:) = ssht_forward(reshape(fg1(temp+1:temp+L*(2*L-1)),L,2*L-1),L);
end
t_shfor = toc;
display('max flm in highest frequency, vol1',num2str(max(max(abs(flm1(:,(L-1)^2+1:end))))));
tic
[ coeff ]= FBcoeff_sphere(flm1, basis, sample_points);
t_fbfor = toc;

% 0. generate random SB coefficients
maxL = L-1; r_mask = R; r_select_ratio = 1;
load SphericalBessel.mat
B = bessel;
B = B( B(:, 4)<= pi * r_mask * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
l_grid = B(:, 1);
n_grid = B(:, 2);
clear B bessel;
siz_n = zeros(maxL+1,1);
for ll = 0: maxL
    ind_ls = find( l_grid == ll);
    siz_n(ll+1) = numel( ind_ls );
end
A = coeff2cell(coeff);

% 3. compute its steerable PCA to get eigenvolumes
Avec = cell2vec(A);
C = Avec*Avec';
[U,S,~] = svd(C,'econ');
eigVec = U(:,1);
eigVal = S(1,1);
clear U S V C;
cellSteer = vec2cell(eigVec, siz_n);
coeffSteer = cell2coeff(cellSteer);


Cl = cell(maxL+1,1);
SVD = cell(maxL+1,3);
eigVecSteer = cell(maxL+1,1);
eigValSteer = cell(maxL+1,1);
for ll = 0:maxL
    Cl{ll+1} = A{ll+1}*A{ll+1}';
    [SVD{ll+1,1}, SVD{ll+1,2}, SVD{ll+1,3}] = svd(Cl{ll+1},'econ');
    % steerable PCA
    ind_max = min(2*ll+1,siz_n(ll+1));
    eig = SVD{ll+1,1};
    eigVecSteer{ll+1} = eig(:,1:ind_max);
    eig = diag(SVD{ll+1,2});
    eigValSteer{ll+1} = eig(1:ind_max);   
end
eigL = cell2vec(eigVecSteer);

%% 4. compute N rotated version of it from coefficients
E = zeros(4,1);
for j = 1:4
N = 10^j;
angles = rand(3,N)*2*pi;
angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
hist(angles(2,:));
A_N = cell(maxL+1, N);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        [T,Tinv] = realY_to_complexY(ll);
        Wl = real(Tinv*Wl*T);
        al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
        A_N{ll+1,i} = al;
    end
end

% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
ANvec = zeros(size(Avec,1),N);
for i = 1:N
    ANvec(:,i) = cell2vec(A_N(:,i));
end
ANvecMean = mean(ANvec,2);
rel = norm([ANvecMean(1:20)-coeff{1,1}; ANvecMean(21:end)])/norm(ANvecMean);
E(j) = rel;
save('error.mat','E');
end

%%
Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))';
[Us,Ss,Vs] = svd(Cs,'econ');
eigValS = Ss(1,1);
eigVecS = Us(:,1);
cellS = vec2cell(Us(:,1), siz_n);
coeffS = cell2coeff(cellS);

% Cls = cell(maxL+1,1);
% SVDs = cell(maxL+1,3);
% eigVecSP = cell(maxL+1,1); % sample with classical PCA
% eigValSP = cell(maxL+1,1);
% for ll = 0:maxL
%     Cls{ll+1,1} = cellfun(@(x) x*x', A_N(ll+1,:), 'UniformOutput',false);%A{ll+1}*A{ll+1}'; curly bracket for content, paren for cell
%     Cls{ll+1,2} = mean(cat(3,Cls{ll+1,1}{:}),3);
%     [SVDs{ll+1,1}, SVDs{ll+1,2}, SVDs{ll+1,3}] = svd(Cls{ll+1,2},'econ');
%     % steerable PCA
%     ind_max = min(2*ll+1,siz_n(ll+1));
%     eig = SVDs{ll+1,1};
%     eigVecSP{ll+1} = eig(:,1:ind_max);
%     eig = diag(SVDs{ll+1,2});
%     eigValSP{ll+1} = eig(1:ind_max);
% end

figure()
subplot(1,2,1)
plot(1:length(ANvecMean),ANvecMean)
legend('Sample mean')
xlabel('coefficient a_{nlm}');
subplot(1,2,2)
plot(1:20,ANvecMean(1:20),'o'); hold on
plot(1:20,coeff{1,1});
xlabel('coefficient a_{nlm}');
legend('Sample mean','Steerable PCA mean');
title('compare sample mean volume and that of steerable PCA as coefficient');

rel = norm([ANvecMean(1:20)-coeff{1,1}; ANvecMean(21:end)])/norm(ANvecMean);

[ vol_hatSMean, err, t_ift, ball ] = IFT_SB(R, c, L, cell2coeff(vec2cell(ANvecMean,siz_n)), vol);
coeffSteerM = zeros(size(ANvecMean));
coeffSteerM(1:length(coeff{1,1})) = coeff{1,1};
[ vol_hatSteerM, err, t_ift, ball ] = IFT_SB(R, c, L, cell2coeff(vec2cell(coeffSteerM,siz_n)), vol);

figure()
subplot(1,2,1)
threshF = quantile(vol_hatSMean(:),0.8);
Fv = isosurface(X,Y,Z,vol_hatSMean,threshF);
p = patch(Fv);
%isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('recon of sample mean in real space');

subplot(1,2,2)
threshF = quantile(vol_hatSteerM(:),0.8);
Fv = isosurface(X,Y,Z,vol_hatSteerM,threshF);
p = patch(Fv);
%isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('recon of steerable PCA mean in real space');


figure()
plot(real(eigVec)); hold on; plot(real(eigVecS)); hold on; plot(real(Us(:,1)));
legend('steerable PCA','1st component of sample','2nd component');

%%
[ vol_hatSteer, err, t_ift, ball ] = IFT_SB(R, c, L, coeffSteer, vol);
figure()
subplot(2,5,1)
threshF = quantile(vol_hatSteer(:),0.1);
Fv = isosurface(X,Y,Z,vol_hatSteer,threshF);
p = patch(Fv);
%isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
title('recon of steerable PCA in real space');

for i = 1:9
    cellS = vec2cell(Us(:,i), siz_n);
    coeffS = cell2coeff(cellS);
    [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, c, L, coeffS, vol);
    subplot(2,5,1+i)
    threshF = quantile(vol_hat1(:),0.1);
    Fv = isosurface(X,Y,Z,vol_hat1,threshF);
    p = patch(Fv);
    %isonormals(xq,yq,zq,vq,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    title([num2str(i) '-th eigenvolume lambda=' num2str(sqrt(Ss(i,i)),3)]);
end


figure()
plot(diag(Ss));
%% 6. compare the eigenvalues and eigenvectors, difference to machine
% precision
% dEVal = cellfun(@(x,y) norm(x-y)/norm(x), eigValSteer, eigValSP, 'UniformOutput',false);
% dEVec = cellfun(@(x,y) norm(x-y)/norm(x), eigVecSteer, eigVecSP, 'UniformOutput',false);

% 7. compute the volume for the eigenvectors

[ dataS ]= FBcoeff_sphere_inv(coeffS, basis, sample_points);
[ dataSteer ] = FBcoeff_sphere_inv(coeffSteer, basis, sample_points);
fS = [];
fSteer = [];
for i = 1:n_r
    tempS = ssht_inverse(dataS(i,:),L);
    tempSteer = ssht_inverse(dataSteer(i,:),L);
    fS = cat(1,fS,tempS(:));
    fSteer = cat(1,fSteer, tempSteer(:));
end


% [thetas,phis] = ssht_sampling(L,'Grid',true);
% [x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));% thetas(:)-pi/2
% omega = [x y z];
% omega_all = omega*sample_points.r(1);
% for i = 2:n_r
%     omega_all = cat(1,omega_all,omega*sample_points.r(i));
% end;

% figure(1);
% subplot(1,2,1)
% hist(real(fS));
% subplot(1,2,2)
% hist(real(fSteer))

FSteer = scatteredInterpolant(omega_all(:,1),omega_all(:,2),omega_all(:,3),real(fSteer));
[xq,yq,zq] = meshgrid(-pi:0.25:pi);
vq = FSteer(xq,yq,zq);
threshSteer = quantile(real(fSteer),0.3);
fv = isosurface(xq,yq,zq,vq,threshSteer);
p = patch(fv);
isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

FS = scatteredInterpolant(omega_all(:,1),omega_all(:,2),omega_all(:,3),real(fS));
vq = FS(xq,yq,zq);
threshS = quantile(real(fS),0.1);
fv = isosurface(xq,yq,zq,vq,threshS);
p = patch(fv);
isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'blue';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

title('reconstructed volumes, red-steerable PCA, blue-top component of sample covariance')
