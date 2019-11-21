addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');

% 1. generate random SB coefficients
maxL = 16; R = maxL; r_mask = R; r_select_ratio = 1;
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

% 2. compute its steerable PCA to get eigenvolumes
%A{1,1}= A{1,1}- A{1,1};
Avec = cell2vec(A);
C = zeros(size(Avec,1));
C_cell = cell(maxL+1,1);
C_cell{1,1} = zeros(siz_n(1));
count = siz_n(1);
for ll = 1:maxL
    n = siz_n(ll+1);
    for m = 1:2*ll+1
        C(count+1:count+n,count+1:count+n) = A{ll+1,1}*A{ll+1,1}'/(2*ll+1);
        count = count+n;
    end
    C_cell{ll+1,1} = A{ll+1,1}*A{ll+1,1}'/(2*ll+1);
end

%% 3. compute sample covariance matrix
l = size(Avec,1);
num = 5;
rep = 1;
E = zeros(num,rep,length(C(:)));
T = zeros(num,rep);
Cs_all = cell(maxL+1,num);
for j = 1:num
    display(j)
    for k = 1:rep
        N = 10^j;
        tic
        angles = rand(3,N)*2*pi;
        angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
        AN = cell(maxL+1, 1);
        ANvec = zeros(size(Avec,1),N);
        for i = 1:N
            angle = angles(:,i);
            for ll = 0:maxL
                Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
                %[T,Tinv] = realY_to_complexY(ll);
                %Wl = real(Tinv*Wl*T);
                al = A{ll+1}*Wl;%conj(Wl); % give complex numbers, use real Y?
                AN{ll+1,1} = al;
            end
            ANvec(:,i) = cell2vec(AN);
        end
        toc

        % compute the sample covariance matrix 
        ANvecMean = mean(ANvec,2);
        Cs = zeros(size(Avec,1));
        tic
        for i = 1:floor(N/l)
            Cs = Cs + (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))*...
                (ANvec(:,(i-1)*l+1:i*l)-repmat(ANvecMean,1,l))';
        end
        if floor(N/l) == 0
            Cs = (ANvec - repmat(ANvecMean,1,N))*(ANvec - repmat(ANvecMean,1,N))';
        else
            Cs = Cs + (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N-i*l))*...
                    (ANvec(:,i*l+1:end)-repmat(ANvecMean,1,N-i*l))';
        end
        
        Cs_cell = cell(maxL+1,1);
        count = 0;
        for ll = 0:maxL
            n = siz_n(ll+1);
            cel = zeros(n);
            for m = 1:2*ll+1
                cel = cel + C(count+1:count+n,count+1:count+n);
                count = count+n;
            end
            Cs_cell{ll+1,1} = real(cel/(2*ll+1));
        end
        %T(j,k) = toc;
        rel = norm(real(Cs)/N-C)/norm(C);
        E(j,k) = rel;
        Cs_all(:,j) = Cs_cell;
        display(rel);
        %save('coeff_convB8_num25rep5.mat','E','T');
    end
end


%%
%load('coeff_convB8_num25rep5.mat','E','T');
for j = 1:num
    display(norm(real(squeeze(E(j,:,:))))/norm(C));
end
    

figure()
subplot(1,2,1)
imagesc(real(Cs(1:80,1:80)/N));
colorbar
subplot(1,2,2)
imagesc(real(C(1:80,1:80)));
colorbar

%%

figure()
load('coeff_convB8_num25rep5.mat','E','T');
set(gca,'FontSize',14)
subplot(1,2,1)
loglog(2.^(1:num),E(1:num),'o-');
title('Relative error of the sample covariance vs steerable covariance');
xlabel('Number of samples');
subplot(1,2,2)
loglog(2.^(1:num),T);
title('Time for computing the sample covariance matrix');
xlabel('Number of samples');

