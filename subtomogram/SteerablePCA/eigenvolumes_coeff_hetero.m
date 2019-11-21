addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');

%% 1. generate random SB coefficients for M different structures
maxL = 8; R = maxL; r_mask = R; r_select_ratio = 1;
load SphericalBessel.mat
B = bessel;
B = B( B(:, 4)<= pi * r_mask * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
l_grid = B(:, 1);
n_grid = B(:, 2);
clear B bessel;
siz_n = zeros(maxL+1,1); % size of Bessel zeros
for ll = 0: maxL
    ind_ls = find( l_grid == ll);
    siz_n(ll+1) = numel( ind_ls );
end
M = 3; % number of heterogeneous structures
A = cell(maxL+1,M);
for ll = 0:maxL
    for m = 1:M
        A{ll+1,m} = 2*rand(siz_n(ll+1),2*ll+1)-1;
        if mod(ll,2) == 1
            A{ll+1,m} = A{ll+1,m}*1i;
        end
    end
end

%% 2. compute its steerable PCA to get eigenvolumes
for m = 1:M % subtract the average as the 0th coefficient? or the combination?
    A{1,m}= A{1,m}- A{1,m};
end
Avec = [];
for m = 1:M
    Avec = [Avec cell2vec(A(:,m))];
end
C = zeros(size(Avec,1));
count = 0;
for ll = 0:maxL
    n = siz_n(ll+1);
    for mm = 1:2*ll+1
        for m = 1:M
            C(count+1:count+n,count+1:count+n) = C(count+1:count+n,count+1:count+n)+...
                A{ll+1,m}*A{ll+1,m}'/(2*ll+1);
        end
        count = count+n;
    end
end

%% 3. compute sample covariance matrix
l = size(Avec,1);
num = 25;
rep = 5;
E = zeros(num,rep);
T = zeros(num,rep);
for j = 1:num
    display(j)
    N = 2^j;
    for k = 1:rep
        AN = cell(maxL+1, 1);
        ANvec = zeros(size(Avec,1),N*M);
        for m = 1:M
            angles = rand(3,N)*2*pi;
            angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
            for i = 1:N
                angle = angles(:,i);
                for ll = 0:maxL
                    Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
                    %[T,Tinv] = realY_to_complexY(ll);
                    %Wl = real(Tinv*Wl*T);
                    al = A{ll+1,m}*conj(Wl); % give complex numbers, use real Y?
                    AN{ll+1,1} = al;
                end
                ANvec(:,i+N*(m-1)) = cell2vec(AN);
            end
        end

        % compute the sample covariance matrix 
        ANvecMean = mean(ANvec,2);
        Cs = zeros(size(Avec,1));
        tic
        for i = 1:floor(N*M/l)
            Cs = Cs + (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))*...
                (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))';
        end
        if floor(N*M/l) == 0
            Cs = (ANvec - repmat(ANvecMean,1,N*M))*(ANvec - repmat(ANvecMean,1,N*M))';
        else
            Cs = Cs + (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N*M-i*l))*...
                    (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N*M-i*l))';
        end

        T(j,k) = toc;
        rel = norm(Cs/(N*M)-C/M)/norm(C/M);
        E(j,k) = rel;
        display(rel);
        save('coeff_convB8_num25rep5M3.mat','E','T','Cs','C');
    end
end

figure()
subplot(1,2,1)
imagesc(abs(Cs(1:25,1:25)/N));
colorbar
subplot(1,2,2)
imagesc(abs(C(1:25,1:25)));
colorbar

%%

figure()
num = 19;
set(gca,'FontSize',14)
subplot(1,2,1)
loglog(2.^(1:num),mean(E(1:num,:),2),'o-');
title('Relative error of the sample covariance vs steerable covariance');
xlabel('log(N), number of samples');
subplot(1,2,2)
loglog(2.^(1:num),mean(T(1:num,:),2),'o-');
title('Time for computing the sample covariance matrix');
xlabel('log(N), number of samples');

%% visualize in real space
N = 2^19;
[U,S,~] = svd(C,'econ');
S = diag(S);
U = U(:,1:100);
% save('eigVols.mat','U','-append');
% save('eigVols.mat','S','-append');

[Us, Ss, ~] = svd(Cs/N,'econ');
Ss = diag(Ss);
Us = Us(:,1:100);
% save('eigVols.mat','Us','-append');
% save('eigVols.mat','Ss','-append');

figure();
subplot(1,2,1)
imagesc(abs(Cs(1:100,1:100))/N);
colorbar
subplot(1,2,2)
imagesc(abs(C(1:100,1:100)));
colorbar

figure(); 
set(gca,'FontSize',14)
loglog(2.^(1:12),E(1:12),'o-');
title('Relative error of the sample covariance vs steerable covariance');
xlabel('log(N), number of samples');

R = 8;
nx = 2*R+1;
x1 = -R:1:R;
F3d = zeros(nx,nx,nx);
[X,Y,Z] = meshgrid(x1,x1,x1);
% Fourier domain
c = pi;
L = R+1; 

%%
figure()
for i = 1:5
    coeffSt = cell2coeff(vec2cell(U(:,i),siz_n));
    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, c, L, coeffSt, F3d);
    subplot(2,5,i)
    threshF = quantile(vol_hatSt(:),0.1);
    Fv = isosurface(X,Y,Z,vol_hatSt,threshF);
    p = patch(Fv);
    %isonormals(xq,yq,zq,vq,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    title(['Steer PCA eigenvol,' num2str(i) num2str(sqrt(S(i)),'%10.3e')]);
end

for i = 1:5
    cellS = vec2cell(Us(:,i), siz_n);
    coeffS = cell2coeff(cellS);
    [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, c, L, coeffS, F3d);
    subplot(2,5,5+i)
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
    title(['PCA eigenvol,' num2str(i) num2str(sqrt(Ss(i)),'%10.3e')]);
end

