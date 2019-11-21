addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/aspire'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ffb'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ssht'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'),...
    genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/easyspin-5.0.20/easyspin'),...
    '/u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan/NUG/realWigner');

load('matlab.mat');
% 3. compute its steerable PCA to get eigenvolumes
A = coeff2cell(coeff);
A{1,1}= A{1,1}- coeff{1,1};
Avec = cell2vec(A);
C = Avec*Avec';

%%

N = 10;
alpha = linspace(0,2*pi,N+1);
gamma = alpha;
beta = linspace(0,pi,N+1);
[x,y,z] = meshgrid(alpha,beta,gamma);
angles = [x(:) y(:) z(:)]';
% angles = rand(3,N)*2*pi;
% angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
hist(angles(2,:));
A_N = cell(maxL+1, (N+1)^3);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        %[T,Tinv] = realY_to_complexY(ll);
        %Wl = real(Tinv*Wl*T);
        al = A{ll+1}*conj(Wl)*sqrt(sin(angle(2))); % give complex numbers, use real Y?
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


figure();
subplot(1,2,1)
imagesc(abs(Cs(1:100,1:100))/N);
subplot(1,2,2)
imagesc(abs(C(1:100,1:100)));


% figure(); 
% set(gca,'FontSize',14)
% loglog(2.^(1:13),E,'o-');
% title('Relative error of the sample mean vs invariant mean');
% xlabel('log(N), number of samples');

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




