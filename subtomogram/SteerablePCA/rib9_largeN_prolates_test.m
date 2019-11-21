addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');
addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/pswf_t_3d_opt/pswf_t_3d');


%% Define expansion parameters
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values correspond for greater oversampling (choose beta=1 for Nyquist sampling, beta=0.5 for x2 oversampling, etc.)
delta = 0.99;   % Truncation parameter (between 0 and 1) - small values improve accuracy, large values use shorter expansions and improve noise spectrum. Approximation error should be small even for values near 1 if the volume is spatially localized. 

% For all parameters - see more details in the function pswf_t_f_3d. 


%% downsample ribosome
load('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/projections/simulation/maps/cleanrib.mat')
R = 9;
L = R;
x1 = -R:R;
nx = length(x1);
vol = real(cryo_downsample(volref,[nx nx nx]));

%% Define test volume parameters
gridSize = nx;  % number of voxels in each dimenison
% t = 0.1;        % Gaussian function standard deviation
nVols = 1;      % Number of volume functions (to test timing, all volumes are equal)

%% Duplicate volumes for timing test
vol = repmat(vol,1,1,1,nVols);

%% Obtain expansion coefficients
tic
coeff = pswf_t_f_3d(vol, beta, delta);          % Obtain coefficients of 3D PSWF expansion - CPU version
% coeff = pswf_t_f_3d_gpu(vol, beta, delta);          % Obtain coefficients of 3D PSWF expansion - GPU version
t_f = toc;
disp(['Computing PSWF coefficients took ',num2str(t_f),' seconds.']);

%% change expansion coefficients to full -m to m for rotation
coeff_orig = coeff;
coeff = coeff{1,1};
maxL = size(coeff,2);
for l = 0:maxL-1
    temp_l = coeff{l+1};
    temp_c = temp_l(:,1);
    for m = 1:l
        col = temp_l(:,m+1);
        temp_c = [real(col)+imag(col)*(-1)^m*1i temp_c col];
    end
    coeff{l+1} = temp_c;
end

siz_n = zeros(maxL,1);
for l = 0:maxL-1
    siz_n(l+1) = size(coeff{l+1},1);
end

size_p = siz_n;
coeff_p = coeff;
maxL_p = maxL;
save('rib9_prolates_basics_Nov4.mat','size_p','coeff_p','maxL_p');
%%
% 3. compute its steerable PCA to get eigenvolumes
A = coeff;%coeff2cell(coeff);
A{1,1}= A{1,1}- coeff{1,1};
Avec = cell2vec(A);
C = zeros(size(Avec,1));
count = 0;
for ll = 0:maxL-1
    n = siz_n(ll+1);
    for m = 1:2*ll+1
        C(count+1:count+n,count+1:count+n) = A{ll+1}*A{ll+1}'/(2*ll+1);
        count = count+n;
    end
end


%% steerable PCA

j = 21;
l2 = 10^3;
    display(j)
    N = 2^j;
    angles = rand(3,N)*2*pi;
    angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
    %hist(angles(2,:));
    AN = cell(maxL+1, 1);
    ANvec = zeros(size(C,1),N);
    tic
    for i = 1:N
        angle = angles(:,i);
        for ll = 0:maxL-1
            Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
            %[T,Tinv] = realY_to_complexY(ll);
            %Wl = real(Tinv*Wl*T);
            al = coeff{ll+1}*conj(Wl); % give complex numbers, use real Y?
            AN{ll+1,1} = al;
        end
        ANvec(:,i) = cell2vec(AN);
    end

    % 5. compute classical PCA with the sample covariance matrix 
    %       and obtain eigenvectors 
    ANvecMean = mean(ANvec,2);
    Cs = zeros(size(ANvec,1));

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
    rel = norm(Cs/N-C)/norm(C);
    display(rel);
    Cs_prolates = Cs/N;
    Tp = toc;
    save('rib2019_G21_prolates_Nov4.mat','Cs_prolates','Tp');

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

%%
% figure(1)
% load('coeff_convB8_num25rep5.mat','E','T');
% load('error2019_G25_prolates.mat');
% load('error2019_G25_bessel.mat');
% set(gca,'FontSize',16)
% 
% num = 25;
% cur = 16;
% E = mean(E,2);
% T = mean(T,2);
% subplot(1,2,1)
% loglog(2.^(1:cur),E2(1:cur),'o-',2.^(1:num),E3(1:num),'o-');
% title('Relative error of sample covariance vs steerable covariance','FontSize',11);
% xlabel('Number of samples');
% ylabel('relative error in 2-norm');
% legend('Spherical Bessel','Prolates');
% subplot(1,2,2)
% loglog(2.^(1:cur),T2(1:cur),'o-',2.^(1:num),T3,'o-');
% title('Time to compute sample covariance matrix for volume of size 25','FontSize',11);
% xlabel('Number of samples');
% ylabel('Time (seconds)');
% legend('Spherical Bessel','Prolates');
% h=gcf;
% set(h,'PaperOrientation','landscape');
% print(gcf,'-fillpage','cov_Gaussian','-dpdf')