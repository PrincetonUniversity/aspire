% 0. generate random SB coefficients
% 3. compute its steerable PCA to get eigenvolumes
% 4. compute N rotated version of it
% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
% 6. compare eigenvectors and eigenvalues from classical and steerable PCA

% Nov 27, 2017
addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner');


% 0. generate random SB coefficients
maxL = 8; R = maxL; r_mask = R; r_select_ratio = 1;
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
A = cell(maxL+1,1);
for ll = 0:maxL
    A{ll+1} = 2*rand(siz_n(ll+1),2*ll+1)-1;
    if mod(ll,2) == 1
        A{ll+1} = A{ll+1}*1i;
    end
end

% 3. compute its steerable PCA to get eigenvolumes
Avec = cell2vec(A);
C = Avec*Avec';
[U,S,V] = svd(C,'econ');

% Cl = cell(maxL+1,1);
% SVD = cell(maxL+1,3);
% eigVecSteer = cell(maxL+1,1);
% eigValSteer = cell(maxL+1,1);
% for ll = 0:maxL
%     Cl{ll+1} = A{ll+1}*A{ll+1}';
%     [SVD{ll+1,1}, SVD{ll+1,2}, SVD{ll+1,3}] = svd(Cl{ll+1},'econ');
%     % steerable PCA
%     ind_max = min(2*ll+1,siz_n(ll+1));
%     eig = SVD{ll+1,1};
%     eigVecSteer{ll+1} = eig(:,1:ind_max);
%     eig = diag(SVD{ll+1,2});
%     eigValSteer{ll+1} = eig(1:ind_max);   
% end
% eigL = cell2vec(eigVecSteer);

%% 4. compute N rotated version of it from coefficients
N = 10000;
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
Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))'; %-repmat(ANvecMean,1,N)
[Us,Ss,Vs] = svd(Cs,'econ');
figure();
plot(diag(Ss));


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
%hold on; plot(-real(Us(1:10,1))); 
%% 
figure()
subplot(2,3,1);
plot(real(Avec(:)),'o'); hold on; 
plot(imag(Avec(:)),'o'); hold on;
plot(sqrt(S(1,1))*real(U(:,1)),'+'); hold on; 
plot(sqrt(S(1,1))*imag(U(:,1)),'+'); hold on; 
legend('original coeff real','original coeff complex',...
    'steerable PCA real','steerable PCA complex');
title('steerable PCA approximation')
subplot(2,3,2);
temp = siz_n(1);
plot(sqrt(S(1,1))*real(U(:,1)),'o-'); hold on; 
fac = Ss(1,1)/sum(diag(Ss(1:maxL+1,1:maxL+1)))*sqrt(S(1,1));
plot(real(Us(:,1)),'x-');
legend('steerable PCA','1st PC');
subplot(2,3,3);
temp_new = siz_n(2);
plot(-sqrt(S(1,1))*imag(U(temp+1:temp+temp_new*3,1)),'o-'); hold on; 
%fac = sqrt(S(1,1))*Ss(1,1)/sum(diag(Ss(1:5,1:5)));
plot(-imag(Us(temp+1:temp+temp_new*3,2)),'x-');
legend('steerable PCA','2nd PC');

figure();
plot(real(Us(:,1))); hold on;
plot(real(U(:,1)));


%% 6. compare the eigenvalues and eigenvectors, difference to machine
% precision
dEVal = cellfun(@(x,y) norm(x-y)/norm(x), eigValSteer, eigValSP, 'UniformOutput',false);
dEVec = cellfun(@(x,y) norm(x-y)/norm(x), eigVecSteer, eigVecSP, 'UniformOutput',false);

% 7. compute the volume for the eigenvectors
cellS = vec2cell(Us(:,1), siz_n);
coeffS = cell2coeff(cellS);
cellSteer = vec2cell(U(:,1), siz_n);
coeffSteer = cell2coeff(cellSteer);
c = pi;
R = maxL;
n_r = ceil(2*c*R); %2cR according to Jane, 4cR Nyquist? or sqrt(3)*c according to SB
[ basis, sample_points ] = precomp_fb_sphere( n_r, R, c, 2, maxL+1); 

[ dataS ]= FBcoeff_sphere_inv(coeffS, basis, sample_points);
[ dataSteer ] = FBcoeff_sphere_inv(coeffSteer, basis, sample_points);
L = maxL+1;
fS = [];
fSteer = [];
for i = 1:n_r
    tempS = ssht_inverse(dataS(i,:),L);
    tempSteer = ssht_inverse(dataSteer(i,:),L);
    fS = cat(1,fS, tempS(:));
    fSteer = cat(1,fSteer, tempSteer(:));
end

[thetas,phis] = ssht_sampling(L,'Grid',true);
[x, y, z] = sph2cart(phis(:),thetas(:)-pi/2,ones(L*(2*L-1),1));% thetas(:)-pi/2
omega = [x y z];
omega_all = omega*sample_points.r(1);
for i = 2:n_r
    omega_all = cat(1,omega_all,omega*sample_points.r(i));
end;

figure(1);
subplot(1,2,1)
hist(real(fS));
subplot(1,2,2)
hist(real(fSteer))
threshS = quantile(real(fS),0.92);
threshSteer = quantile(real(fSteer),0.92);

FS = scatteredInterpolant(omega_all(:,1),omega_all(:,2),omega_all(:,3),real(fSteer));
[xq,yq,zq] = meshgrid(-pi:0.25:pi);
vq = FS(xq,yq,zq);
fv = isosurface(xq,yq,zq,vq,threshS);
p = patch(fv);
isonormals(xq,yq,zq,vq,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud





[thetas,phis] = ssht_sampling(L,'Grid',true);
type = 'colour';
figure(1);
subplot(1,2,1)
ssht_plot_sphere(real(tempS), L);
title('top eigenvector from sample');
subplot(1,2,2)
ssht_plot_sphere(real(tempSteer),L);
title('Steerable PCA eigenvector');