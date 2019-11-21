addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');

load('matlab.mat');

%%
% 3. compute its steerable PCA to get eigenvolumes
A = coeff2cell(coeff);
A{1,1}= A{1,1}- coeff{1,1};
C = zeros(size(cell2vec(A),1));
count = 0;
for ll = 0:maxL
    n = siz_n(ll+1);
    for m = 1:2*ll+1
        C(count+1:count+n,count+1:count+n) = A{ll+1,1}*A{ll+1,1}'/(2*ll+1);
        count = count+n;
    end
end

num = 13;
E = zeros(num,1);
for j = 12:12
display(j)
N = 2^j;
angles = rand(3,N)*2*pi;
angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
hist(angles(2,:));
A_N = cell(maxL+1, N);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        %[T,Tinv] = realY_to_complexY(ll);
        %Wl = real(Tinv*Wl*T);
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
Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))';
rel = norm(Cs/N-C)/norm(C);
E(j) = rel;
display(rel);
%save('errorC2.mat','E');
end

save('eigVols.mat','C','Cs');
[U,S,~] = svd(C,'econ');
S = diag(S);
U = U(:,1:100);
save('eigVols.mat','U','-append');
save('eigVols.mat','S','-append');

[Us, Ss, ~] = svd(Cs/N,'econ');
Ss = diag(Ss);
Us = Us(:,1:100);
save('eigVols.mat','Us','-append');
save('eigVols.mat','Ss','-append');

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


%% confirm the orthogonality of wignerD matrix
E = zeros(9,1);
for j = 2:10
N = 2^j;
angles = rand(3,N)*2*pi;
angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
hist(angles(2,:),20);
WL = cell(maxL+1,N);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:1
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        WL{ll+1,i} = Wl;
    end
end
% integrate conj(W^l_{mm})*W^l_{mm}
l1 = 1;
l2 = 1;
m11 = 3;
m12 = 3;
m21 = 2;
m22 = 3;
int1 = cellfun(@(x1,x2) conj(x1(m11,m12))*x2(m21,m22), WL(l1+1,:), ...
    WL(l2+1,:), 'UniformOutput',false);
E(j-1) = sum(cell2mat(int1).*sin(angles(2,:)))/N;
end




