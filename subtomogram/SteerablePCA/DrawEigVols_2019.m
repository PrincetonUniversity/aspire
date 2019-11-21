% align the principal components, invert them back to the real domain,
% and draw central slices

% EigVols with 2^15 samples
% load('eigVols.mat');
% load('rib_basics.mat');
% i = 1;
% vecSt = U(:,i);
% vecS = Us(:,i);
% %eul1 = alignComp(vecSt,vecS,siz_n);
% % eul1(2) = -(180-eul1(2));
% cell1 = vec2cell(vecSt,siz_n);
% maxL = size(cell1,1)-1;
% Ctemp = cell(1,maxL+1);
% Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
% Wp = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
% Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
% Wm = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
% cellSR = rotateCellFast(cell1,[pi/4,0,0]',Wp,Wm);
% eul1 = alignComp(cell2vec(cellSR),vecSt, siz_n);

%%
% L = 16;
% N = 2*L-1;
% alpha = pi*(0:N)/L;
% beta = pi*((0:N)*2+1)/(4*L);
% [alpha,beta,gamma] = ndgrid(alpha,beta,alpha);
% angles = [alpha(:) beta(:) gamma(:)];
% A = vec2cell(vecSt,siz_n);
% AN = cell(maxL+1, 1);
% ANvec = zeros(size(Avec,1),(2*L)^3);
% for i = 1:(2*L)^3
%     angle = angles(i,:)';
%     for ll = 0:maxL
%         Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
%         %[T,Tinv] = realY_to_complexY(ll);
%         %Wl = real(Tinv*Wl*T);
%         al = A{ll+1,1}*conj(Wl); % give complex numbers, use real Y?
%         AN{ll+1,1} = al;
%     end
%     ANvec(:,i) = cell2vec(AN);
% end
% diff = zeros((2*L)^3,1);
% for i = 1:(2*L)^3
%     diff(i) = norm(ANvec(:,i)-vecS);
% end
% [~,ind] = min(diff);
% 
% 
% %%
% num = 5;
% L = 9;
% figure()
% for i = 1%:num
%     vecS = Us(:,i);
% %     vecS = Us(:,i);
% %     eul1 = alignComp(vecS,vecSt,siz_n);
% %     eul1(2) = -eul1(2);
% %     cellSR = rotateCell(vec2cell(vecSt,siz_n),eul1./180);
% %     coeffSt = cell2coeff(cellSR);
%     coeffSt = vec2cell(vecS,siz_n);%cell2coeff
%     
%     [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, c, L, coeffSt, F3d);
%     subplot(2,num,i)
%     threshF = quantile(vol_hatSt(:),0.1);
%     Fv = isosurface(X,Y,Z,vol_hatSt,threshF);
%     p = patch(Fv);
%     %isonormals(xq,yq,zq,vq,p)
%     p.FaceColor = 'red';
%     p.EdgeColor = 'none';
%     daspect([1 1 1])
%     view(3); 
%     axis tight
%     camlight 
%     lighting gouraud
%     title(['Steer PCA eigenvol,' num2str(i) num2str(sqrt(S(i)),'%10.3e')]);
% end
% 
% for i = 1%:num
% %     cellS = vec2cell(Us(:,i), siz_n);
% %     coeffS = cell2coeff(cellS);
%     coeffS = cell2coeff(vec2cell(ANvec(:,ind),siz_n));
%     [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, c, L, coeffS, F3d);
%     subplot(2,num,num+i)
%     threshF = quantile(vol_hat1(:),0.1);
%     Fv = isosurface(X,Y,Z,vol_hat1,threshF);
%     p = patch(Fv);
%     %isonormals(xq,yq,zq,vq,p)
%     p.FaceColor = 'red';
%     p.EdgeColor = 'none';
%     daspect([1 1 1])
%     view(3); 
%     axis tight
%     camlight 
%     lighting gouraud
%     title(['PCA eigenvol,' num2str(i) num2str(sqrt(Ss(i)),'%10.3e')]);
% end

%% compute block-wise PCA
% [Us, Ss, ~] = svd(Cs/N,'econ');
% Ss = diag(Ss);
% Us = Us(:,1:100);
% save('eigVols.mat','Us','-append');
% save('eigVols.mat','Ss','-append');

% figure();
% subplot(1,2,1)
% imagesc(abs(Cs(1:100,1:100))/N);
% colorbar
% subplot(1,2,2)
% imagesc(abs(C(1:100,1:100)));
% colorbar

load('r9_bessel_noS.mat');
load('rib9_bessel_basics_noS.mat');
% load('eigVols.mat');
% load('rib_basics.mat');

L = maxL+1;
cov_steer = cell(L,3);
count = 0;

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
    l = Ss_ln(i,1)-1;
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
load('rib9_cPCA_Nov.mat','Cs_pca');
[Uc, Sc, ~] = svd(Cs_pca,'econ');
save('rib9_cPCA_Nov_decomp.mat','Uc','Sc');


%%
%load('eigVols.mat');
Ctemp = cell(1,L);
Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
Wp = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);
Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
Wm = cellfun(@wignerd, num2cell(0:L-1), Ctemp, 'UniformOutput', false);

num = 5; R = 9; c = pi;
mid = R+1;
figure()
% bessel steerable PCA
for i = 1:num
%     vecSt = vecS(:,i);
%     vecSs = Us(:,i);
%     eul1 = alignComp(vecSs,vecSt,siz_n);
%     eul1(2) = -eul1(2);
%     cellSR = rotateCellFast(vec2cell(vecSt,siz_n),eul1'./180,Wp,Wm);
    %coeffSt = cell2coeff(cellSR);
    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, c, L, ...
        vec2cell(vecS(:,i),siz_n));
    subplot(3,num,i)
    imagesc(vol_hatSt(:,:,mid)); colormap(flipud(gray)); axis square; 
    %colorbar;
    title(['Bessel eigenval:' num2str(sqrt(Ss_top(i)),'%10.1e'), ...
        ',l=',num2str(Ss_ln(i,1))],'FontSize',8);
end
num = 5; R = 9; c = 1/2;
mid = R+1;
% prolates steerable PCA
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values correspond for greater oversampling (choose beta=1 for Nyquist sampling, beta=0.5 for x2 oversampling, etc.)
delta = 0.99;   % Truncation parameter (between 0 and 1) - small values improve accuracy, large values use shorter expansions and improve noise spectrum. Approximation error should be small even for values near 1 if the volume is spatially localized. 
gridSize = L;
for i = 1:num
    cellp = transpose(vec2cell(vecSp(:,i),size_p));
    for l = 1:L
        temp = cellp{l};
        cellp{l} = temp(:,l:end);
    end
    vol_hat = pswf_t_b_3d({cellp}, gridSize, beta, delta);       % Reconstruct volume function from 3D PSWF expansion coefficients - CPU version
    
    subplot(3,num,num+i)
    imagesc(vol_hat(:,:,mid)); colormap(flipud(gray)); axis square;
    %colorbar;
    title(['Prolates eigenval:' num2str(sqrt(Sp_top(i)/2^21),'%10.1e'), ...
        ',l=',num2str(Sp_ln(i,1))],'FontSize',8);
end

% classical PCA
num = 5;
nx = 2*9+1;
mid = 10;
for i = 1:num
    vol_hat1 = reshape(Uc(:,i), [nx,nx,nx]);
    subplot(3,num,2*num+i)
    imagesc(vol_hat1(:,:,mid)); colormap(flipud(gray)); axis square;
    title(['PCA eigenval:'  num2str(sqrt(Sc(i)),'%10.1e')],'FontSize',8);
end
h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'-fillpage','eigenvol_rib9_Nov','-dpdf')
