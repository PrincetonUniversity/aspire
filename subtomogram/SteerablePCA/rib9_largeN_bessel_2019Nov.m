addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner',...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/slepian_alpha-master');
addpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/pswf_t_3d_opt/pswf_t_3d');


%% downsample ribosome
load('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/projections/simulation/maps/cleanrib.mat')
R = 9;
L = R;
x1 = -R:R;
nx = length(x1);
vol = real(cryo_downsample(volref,[nx nx nx]));

%% Obtain expansion coefficients
c = 1/2;
ss = 2;
n_r = 2*c*R;
% precompute
%[ basis, sample_points, L ] = precomp_fb_sphere(n_r, R, c, ss); 
[ coeff, basis, sample_points, L] = spherical_expansion_ss(vol, n_r, ss);

%% 0. generate SB coefficients
maxL = L-1; r_mask = R; r_select_ratio = 1;
load SphericalBesselL700R300pi.mat
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

%% change expansion coefficients to full -m to m for rotation
% siz_n = zeros(maxL,1);
% for l = 0:maxL-1
%     siz_n(l+1) = size(coeff{l+1},1);
% end

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

save('rib9_bessel_basics_Nov4.mat','siz_n','coeff','maxL','C');

%% compute steerable PCA in the coefficient space

j = 21;
l2 = 10^3;
    disp(j);
    N = 2^j;
    angles = rand(3,N)*2*pi;
    angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
    %hist(angles(2,:));
    AN = cell(maxL+1, 1);
    ANvec = zeros(size(C,1),N);
    tic
    for i = 1:N
        angle = angles(:,i);
        for ll = 0:L-1
            Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
            %[T,Tinv] = realY_to_complexY(ll);
            %Wl = real(Tinv*Wl*T);
            al = coeff{ll+1}*conj(Wl); % give complex numbers, use real Y?
            AN{ll+1,1} = al;
        end
        ANvec(:,i) = cell2vec(AN);
    end

    % 5. compute the sample covariance matrix 
    %       and obtain eigenvectors 
    ANvecMean = mean(ANvec,2);
    Cs = zeros(size(ANvec,1));
    
    for i = 1:floor(N/l2)
        Cs = Cs + (ANvec(:,(i-1)*l2+1:i*l2)-repmat(ANvecMean,1,l2))*...
            (ANvec(:,(i-1)*l2+1:i*l2)-repmat(ANvecMean,1,l2))';
        if mod(i, 100) == 0
            disp(i);
        end
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
    Cs_bessel = Cs/N;
    T = toc;
    save('r9_bessel_Nov4.mat','Cs_bessel','T');

% figure(1)
% load('coeff_convB8_num25rep5.mat','E','T');
% load('error2019_G25_prolates.mat');
% load('error2019_G25_bessel.mat');
% 
% num = 25;
% cur = 14;
% E = mean(E,2);
% T = mean(T,2);
% set(gca,'FontSize',16)
% subplot(1,2,1)
% loglog(2.^(1:cur),E2(1:cur),'o-',2.^(1:num),E3(1:num),'o-');
% title('Relative error of sample covariance vs steerable covariance');
% xlabel('Number of samples');
% ylabel('relative error in 2-norm');
% legend('Spherical Bessel','Prolates');
% subplot(1,2,2)
% loglog(2.^(1:cur),T2(1:cur),'o-',2.^(1:num),T3,'o-');
% title('Time to compute the sample covariance matrix');
% xlabel('Number of samples');
% ylabel('Time (seconds)');
% legend('Spherical Bessel','Prolates');
% h=gcf;
% set(h,'PaperOrientation','landscape');
% print(gcf,'-fillpage','cov_Gaussian','-dpdf')
