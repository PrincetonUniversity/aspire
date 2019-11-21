% 0. generate random SB coefficients
% 3. compute its steerable PCA to get eigenvolumes
% 4. compute N rotated version of it
% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
% 6. compare eigenvectors and eigenvalues from classical and steerable PCA

% Nov 27, 2017


% 0. generate random SB coefficients
maxL = 16; R = 16; r_mask = R; r_select_ratio = 1;
load SphericalBessel.mat
B = bessel;
B = B( B(:, 3)<pi * r_mask * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
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

% 3. compute its steerable PCA to get eigenvolumes
C = cell(maxL+1,1);
SVD = cell(maxL+1,3);
eigVecSteer = cell(maxL+1,1);
eigValSteer = cell(maxL+1,1);
for ll = 0:maxL
    C{ll+1} = A{ll+1}*A{ll+1}';
    [SVD{ll+1,1}, SVD{ll+1,2}, SVD{ll+1,3}] = svd(C{ll+1},'econ');
    % steerable PCA
    ind_max = min(2*ll+1,siz_n(ll+1));
    eig = SVD{ll+1,1};
    eigVecSteer{ll+1} = eig(:,1:ind_max);
    eig = diag(SVD{ll+1,2});
    eigValSteer{ll+1} = eig(1:ind_max);   
end

% 4. compute N rotated version of it from coefficients
N = 10;
angles = rand(3,N)*2*pi;
angles(2,:) = angles(2,:)/2;
A_N = cell(maxL+1, N);
for i = 1:N
    angle = angles(:,i);
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        [T,Tinv] = realY_to_complexY(ll);
        Wl = real(Tinv*Wl*T);
        al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
        A_N{ll+1,i} = al;
    end
end

% 5. compute classical PCA with the sample covariance matrix 
%       and obtain eigenvectors 
Cs = cell(maxL+1,1);
SVDs = cell(maxL+1,3);
eigVecSP = cell(maxL+1,1); % sample with classical PCA
eigValSP = cell(maxL+1,1);
for ll = 0:maxL
    Cs{ll+1,1} = cellfun(@(x) x*x', A_N(ll+1,:), 'UniformOutput',false);%A{ll+1}*A{ll+1}'; curly bracket for content, paren for cell
    Cs{ll+1,2} = mean(cat(3,Cs{ll+1,1}{:}),3);
    [SVDs{ll+1,1}, SVDs{ll+1,2}, SVDs{ll+1,3}] = svd(Cs{ll+1,2},'econ');
    % steerable PCA
    ind_max = min(2*ll+1,siz_n(ll+1));
    eig = SVDs{ll+1,1};
    eigVecSP{ll+1} = eig(:,1:ind_max);
    eig = diag(SVDs{ll+1,2});
    eigValSP{ll+1} = eig(1:ind_max);
end

% 6. compare the eigenvalues and eigenvectors, difference to machine
% precision
dEVal = cellfun(@(x,y) norm(x-y)/norm(x), eigValSteer, eigValSP, 'UniformOutput',false);
dEVec = cellfun(@(x,y) norm(x-y)/norm(x), eigVecSteer, eigVecSP, 'UniformOutput',false);
