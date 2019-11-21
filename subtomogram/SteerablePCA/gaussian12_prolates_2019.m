addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');

% load('rib3.mat');
% load('rib_basics.mat');
load('matlab_Gaussian12.mat');
%%
% 3. compute its steerable PCA to get eigenvolumes
A = coeff;%coeff2cell(coeff);
A{1,1}= A{1,1}- coeff{1,1};
Avec = cell2vec(A);
C = zeros(size(Avec,1));
count = 0;
maxL = L-1;
for ll = 0:maxL
    n = siz_n(ll+1);
    for m = 1:2*ll+1
        C(count+1:count+n,count+1:count+n) = A{ll+1,1}*A{ll+1,1}'/(2*ll+1);
        count = count+n;
    end
end

%l = size(Avec,1);
num = 21;
%% steerable PCA
% E = zeros(num,1);
% T = zeros(num,1);
% for j = 1:num
%     display(j)
%     N = 2^j;
%     angles = rand(3,N)*2*pi;
%     angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
%     %hist(angles(2,:));
%     AN = cell(maxL+1, 1);
%     ANvec = zeros(size(Avec,1),N);
%     for i = 1:N
%         angle = angles(:,i);
%         for ll = 0:maxL
%             Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
%             %[T,Tinv] = realY_to_complexY(ll);
%             %Wl = real(Tinv*Wl*T);
%             al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
%             AN{ll+1,1} = al;
%         end
%         ANvec(:,i) = cell2vec(AN);
%     end
% 
%     % 5. compute classical PCA with the sample covariance matrix 
%     %       and obtain eigenvectors 
%     ANvecMean = mean(ANvec,2);
%     Cs = zeros(size(Avec,1));
%     tic
%     for i = 1:floor(N/l)
%         Cs = Cs + (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))*...
%             (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))';
%     end
%     if floor(N/l) == 0
%         Cs = (ANvec - repmat(ANvecMean,1,N))*(ANvec - repmat(ANvecMean,1,N))';
%     else
%         Cs = Cs + (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N-i*l))*...
%                 (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N-i*l))';
%     end
% 
%     %Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))';
%     T(j) = toc;
%     rel = norm(Cs/N-C)/norm(C);
%     E(j) = rel;
%     display(rel);
%     save('error2019.mat','E','T');
% end

%% compute classical PCA in the real space

num = 25;
E2 = zeros(num,1);
T2 = zeros(num,1);
dreal = (2*R+1)^3;
l2 = 150;
for j = 1:num
    display(j)
    N = 2^j;
    angles = rand(3,N)*2*pi;
    angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
    %hist(angles(2,:));
    AN = cell(maxL+1, 1);
    ANvec = zeros(dreal,N);
    for i = 1:N
        angle = angles(:,i);
        for ll = 0:maxL
            Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
            %[T,Tinv] = realY_to_complexY(ll);
            %Wl = real(Tinv*Wl*T);
            al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
            AN{ll+1,1} = al;
        end
        [vol_hat,~,~,~] = IFT_SB(R, c, L, AN);
        ANvec(:,i) = vol_hat(:);
    end

    % 5. compute classical PCA with the sample covariance matrix 
    %       and obtain eigenvectors 
    ANvecMean = mean(ANvec,2);
    Cs = zeros(size(ANvec,1));
    tic
    for i = 1:floor(N/l2)
        Cs = Cs + (ANvec(:,(i-1)*l2+1:i*l2)-repmat(ANvecMean,1,l2))*...
            (ANvec(:,(i-1)*l2+1:i*l2)-repmat(ANvecMean,1,l2))';
    end
    if floor(N/l2) == 0
        Cs = (ANvec - repmat(ANvecMean,1,N))*(ANvec - repmat(ANvecMean,1,N))';
    else
        Cs = Cs + (ANvec(:,i*l2+1:end)-repmat(ANvecMean,1,N-i*l2))*...
                (ANvec(:,i*l2+1:end)-repmat(ANvecMean,1,N-i*l2))';
    end

    %Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))';
    T2(j) = toc;
%     rel = norm(Cs/N-C)/norm(C);
%     E2(j) = rel;
%     display(rel);
    Cpca = Cs/N;
    save('error2019_G25.mat','T2','Cpca');
end

% save('eigVols.mat','C','Cs');
% [U,S,~] = svd(C,'econ');
% S = diag(S);
% U = U(:,1:100);
% save('eigVols.mat','U','-append');
% save('eigVols.mat','S','-append');
% 
% [Us, Ss, ~] = svd(Cs/N,'econ');
% Ss = diag(Ss);
% Us = Us(:,1:100);
% save('eigVols.mat','Us','-append');
% save('eigVols.mat','Ss','-append');

