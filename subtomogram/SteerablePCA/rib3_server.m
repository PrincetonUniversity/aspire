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

num = 13;
E = zeros(num,1);
for j = 10:10
display(j)
N = 2^j;
angles = rand(3,N)*2*pi;
angles(2,:) = acos(angles(2,:)/pi-1); % [-1,1]
%hist(angles(2,:));
A_N = cell(maxL+1, N);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        [T,Tinv] = realY_to_complexY(ll);
        Wl = real(Tinv*Wl*T);
        al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
        A_N{ll+1,i} = al;
        if ll == 0
            A_N{ll+1,i} = A_N{ll+1,i} - coeff{1,1};
        end
    end
end

% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
Cs = zeros(size(Avec,1));
%ANvecMean = mean(ANvec,2);
%ANvec = zeros(size(Avec,1),N);
for i = 1:N
    ANvec = cell2vec(A_N(:,i));
    Cs = Cs + ANvec*ANvec'/N;
end

%Cs = (ANvec-repmat(ANvecMean,1,N))*(ANvec-repmat(ANvecMean,1,N))';
rel = norm(Cs-C)/norm(C);
E(j) = rel;
display(rel);
save('errorC1.mat','E');
end


% figure(); 
% set(gca,'FontSize',14)
% loglog(2.^(1:13),E,'o-');
% title('Relative error of the sample mean vs invariant mean');
% xlabel('log(N), number of samples');
