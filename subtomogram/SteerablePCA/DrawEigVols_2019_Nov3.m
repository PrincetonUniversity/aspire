% align the principal components, invert them back to the real domain,
% and draw central slices

addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');
addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/pswf_t_3d_opt/pswf_t_3d');

load('r9_bessel_Nov.mat');
load('rib9_bessel_basics_Nov.mat');
% load('eigVols.mat');
% load('rib_basics.mat');

L = maxL+1;
cov_steer = cell(L,3);
count = 0;
%Cs_bessel = C * C;

for ll = 1:L
    Nl = siz_n(ll);
    cov_steer{ll,1} = Cs_bessel(count+1:count+Nl,count+1:count+Nl);
    [cov_steer{ll,2},cov_steer{ll,3},~] = svd(cov_steer{ll,1},'econ');
    count = count + Nl*(2*ll-1);
end

num = 5;
Ss = [];
for ll = 1:L
    Ss = [Ss; diag(cov_steer{ll,3})];
end
[Ss_top, Ss_ind] = maxk(Ss, num);
Ss_ln = zeros(num,2);
Siz_sum = [0;cumsum(siz_n)];
for i = 1:num
    for ll = 2:L
        if Ss_ind(i) < Siz_sum(ll)
            Ss_ln(i,1) = ll-2;
            Ss_ln(i,2) = Ss_ind(i) - Siz_sum(ll-1);
            break;
        end
    end
end

d = size(Cs_bessel,1);
vecS = zeros(d,num);
for i = 1:num
    l = Ss_ln(i,1);%-1;
    n = Ss_ln(i,2);
    count = 0;
    for ll = 0:l
        count = count + siz_n(ll+1)*(2*ll+1);
    end
    Nl = siz_n(l+1);
    temp = cov_steer{l+1,2};
    for m = -l:l
        vecS(count+1:count+Nl,i) = temp(:,n);
        count = count+Nl;
    end
end

%% plot eigenvolumes from prolates
load('rib2019_G21_prolates.mat');
load('rib9_prolates_basics.mat','size_p','coeff_p','maxL_p');
cov_steer_p = cell(maxL_p,3);
count = 0;

for ll = 1:maxL_p
    Nl = size_p(ll);
    cov_steer_p{ll,1} = Cs3(count+1:count+Nl,count+1:count+Nl);
    [cov_steer_p{ll,2},cov_steer_p{ll,3},~] = svd(cov_steer_p{ll,1},'econ');
    count = count + Nl*(2*ll-1);
end

num = 5;
Sp = [];
for ll = 1:maxL_p
    Sp = [Sp; diag(cov_steer_p{ll,3})];
end
[Sp_top, Sp_ind] = maxk(Sp, num);
Sp_ln = zeros(num,2);
Siz_sum = [0;cumsum(size_p)];
for i = 1:num
    for ll = 2:L
        if Sp_ind(i) < Siz_sum(ll)
            Sp_ln(i,1) = ll-2;
            Sp_ln(i,2) = Sp_ind(i) - Siz_sum(ll-1);
            break;
        end
    end
end

d = size(Cs3,1);
vecSp = zeros(d,num);
for i = 1:num
    l = Sp_ln(i,1);
    n = Sp_ln(i,2);
    count = 0;
    for ll = 0:l-1
        count = count + size_p(ll+1)*(2*ll+1);
    end
    Nl = size_p(l+1);
    temp = cov_steer_p{l+1,2};
    for m = -l:l
        vecSp(count+1:count+Nl,i) = temp(:,n);
        count = count+Nl;
    end
end

%% eigenvolumes for classical PCA
%load('rib9_cPCA_Nov.mat','Cs_pca');
%[Uc, Sc, ~] = svd(Cs_pca,'econ');
%save('rib9_cPCA_Nov_decomp.mat','Uc','Sc');
load('rib9_cPCA_Nov_decomp.mat');

%%
%load('eigVols.mat');
Ctemp = cell(1,L);
Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
Wp = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);
Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
Wm = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);

num = 5; R = 9; c = pi;
mid = R+1;
% figure()
% bessel steerable PCA
%for i = 1:num
%%     vecSt = vecS(:,i);
%%     vecSs = Us(:,i);
%%     eul1 = alignComp(vecSs,vecSt,siz_n);
%%     eul1(2) = -eul1(2);
%%     cellSR = rotateCellFast(vec2cell(vecSt,siz_n),eul1'./180,Wp,Wm);
%    %coeffSt = cell2coeff(cellSR);
%    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, pi, L, ...
%        vec2cell(vecS(:,i),siz_n));
%    subplot(3,num,i)
%    imagesc(vol_hatSt(:,:,mid)); colormap(flipud(gray)); axis square; 
%    %colorbar;
%    title(['Bessel eigenval:' num2str(sqrt(Ss_top(i)),'%10.1e'), ...
%        ',l=',num2str(Ss_ln(i,1))],'FontSize',8);
%end
num = 5; R = 9; c = 1/2;
mid = R+1;
ss = 2;
n_r = 2*c*R;
nx = 2*R+1;
Sc_diag = diag(Sc);
% prolates steerable PCA
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values correspond for greater oversampling (choose beta=1 for Nyquist sampling, beta=0.5 for x2 oversampling, etc.)
delta = 0.99;   % Truncation parameter (between 0 and 1) - small values improve accuracy, large values use shorter expansions and improve noise spectrum. Approximation error should be small even for values near 1 if the volume is spatially localized. 
gridSize = L;

figure();
for i = 1:num
    % Bessel
    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, pi, L, vec2cell(vecS(:,i),siz_n));
    subplot(3,num,i)
    imagesc(vol_hatSt(:,:,mid)); colormap(flipud(gray)); axis square; 
    title(['Bessel eigenval:' num2str(sqrt(Ss_top(i)),'%10.1e'), ...
        ',l=',num2str(Ss_ln(i,1))],'FontSize',8);

    % Prolates
    cellp = transpose(vec2cell(vecSp(:,i),size_p));
    for l = 1:L
        temp = cellp{l};
        cellp{l} = temp(:,l:end);
    end
    vol_hat = pswf_t_b_3d({cellp}, gridSize, beta, delta);       % Reconstruct volume function from 3D PSWF expansion coefficients - CPU version
    [ coeff, basis, sample_points, L] = spherical_expansion_ss(vol_hat, n_r, ss);
    vecSt = vecS(:,i);
    vecSs = cell2vec(coeff);
    eul1 = alignComp(vecSt,vecSs,siz_n);
    eul1(2) = -eul1(2);
    cellSR = rotateCellFast(vec2cell(vecSs,siz_n),eul1'./180,Wp,Wm);
    [ vol_hat, err, t_ift, ball ] = IFT_SB(R, pi, L, cellSR);
    subplot(3,num,num+i)
    imagesc(vol_hat(:,:,mid)); colormap(flipud(gray)); axis square;
    title(['Prolates eigenval:' num2str(sqrt(Sp_top(i)/2^21),'%10.1e'), ...
        ',l=',num2str(Sp_ln(i,1))],'FontSize',8);

    % classical PCA
    vol_hat1 = reshape(Uc(:,i), [nx,nx,nx]);
    [ coeff, basis, sample_points, L] = spherical_expansion_ss(vol_hat1, n_r, ss);
    vecSt = vecS(:,i);
    vecSs = cell2vec(coeff);
    eul1 = alignComp(vecSt,vecSs,siz_n);
    eul1(2) = -eul1(2);
    cellSR = rotateCellFast(vec2cell(vecSs,siz_n),eul1'./180,Wp,Wm);
    [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, pi, L, cellSR);
    subplot(3,num,2*num+i)
    imagesc(vol_hat1(:,:,mid)); colormap(flipud(gray)); axis square;
    title(['PCA eigenval:'  num2str(sqrt(Sc_diag(i)),'%10.1e')],'FontSize',8);
    
end
h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'-fillpage','eigenvol_rib9_Nov10','-dpdf')

% classical PCA
%num = 5;
%nx = 2*9+1;
%mid = 10;
%Sc_diag = diag(Sc);
%for i = 1:num
%    vol_hat1 = reshape(Uc(:,i), [nx,nx,nx]);
%    subplot(3,num,2*num+i)
%    imagesc(vol_hat1(:,:,mid)); colormap(flipud(gray)); axis square;
%    title(['PCA eigenval:'  num2str(sqrt(Sc_diag(i)),'%10.1e')],'FontSize',8);
%end
%h=gcf;
%set(h,'PaperOrientation','landscape');
%print(gcf,'-fillpage','eigenvol_rib9_Nov6','-dpdf')
