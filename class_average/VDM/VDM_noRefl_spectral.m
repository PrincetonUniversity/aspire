function  [r_VDM, r_DM, VV1, id_new]=VDM_noRefl_spectral( rot_matrix, idx_matrix, num_eig, method)
%VDM and DM embedding
%   Input:
%       rot_matrix: rotation for pairs of images of size Px1. in degrees.
%       idx_matrix: list for pairs of size Px2.
%       num_eig: the number of eigenvectors to use.
%       method: 0 means union rule. 1 means intersection rule
%   Output:
%       r_VDM: the embedding from VDM
%       r_DM: the embedding from DM
%       VV1: the top eigenvectors of size N x num_eig.
%   Zhizhen Zhao Feb 10 2012

n_theta=360; %angular discretaizations, 360 radial lines.
N=max(idx_matrix(:)); % the number of images
rows=idx_matrix(:, 1);
cols=idx_matrix(:, 2);

if method==0
    % i and j are nearest neighbors if i is j's nearest neighbor or j is
    % i's nearest neighbors
    C=[rows,cols];
    CC=[cols, rows];
    [~,ia, ib]=union(C, CC, 'rows');
    rows=[C(ia, 1); CC(ib, 1)];
    cols=[C(ia, 2); CC(ib, 2)];
    rot_matrix=[rot_matrix(ia), -rot_matrix(ib)];
    
    ind = find(rows>cols);
    rows=rows(ind);
    cols=cols(ind);
    rot_matrix=rot_matrix(ind);
    A=(rot_matrix)*2*pi/n_theta;
    % using the rot angle and also adding weights
    Ah=exp(-sqrt(-1)*(A));
else
    % i and j are nearest neighbors if i is j's nearest neighbor and j is
    % i's nearest neighbor
    C=[rows,cols];
    CC=[cols, rows];
    [~, ia, ~]=intersect(C, CC, 'rows');
    rows_int=rows(ia);
    cols_int=cols(ia);
    rot=rot_matrix(ia);

    ind=find(rows_int<cols_int);
    rows=rows_int(ind);
    cols=cols_int(ind);
    rot_matrix=rot(ind);
    clear ind
    
    A=(rot_matrix)*2*pi/n_theta;
    
    % using the rot angle and also adding weights
    Ah=exp(-sqrt(-1)*(A));
end;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Building the Class Averaging matrix
H_1=sparse(rows, cols, Ah, N, N);

H_b=H_1+H_1';  % the matrix with original images should be hermitian.

H11=sparse(rows, cols, ones(length(rows), 1), N, N);

H=H11+H11';
% H2_1=sparse(rowsA, colsA, AAh, N, N);
% H2_2=sparse(rowsA, colsA, BBh, N, N);
C=sum(H,2); 
figure; hist(C, 50)
min(C) %test if some nodes are disconnected.
id_new=find(C~=0);
M=length(id_new)
H=H(id_new, id_new);
H_b=H_b(id_new, id_new);
%rebuild H matrix to avoid rows with too less entries
C=sum(H,2); 
figure; hist(C)
A=sparse(1:M, 1:M, 1./sqrt(C), M, M);
H_b1=A*H_b*A;
H1=A*H*A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eigenvector calculation
[VV1, D1]=eigs(H_b1, num_eig);
[VV2, D2]=eigs(H1, num_eig+1);
eig_vals1=diag(D1);
eig_vals2=diag(D2);
[sorted_eigs1, sort_idx1] = sort(real(eig_vals1),'descend');
[sorted_eigs2, sort_idx2] = sort(real(eig_vals2),'descend');
figure; bar(1-sorted_eigs1);
figure; bar(1-sorted_eigs2);
V4=VV1(:, sort_idx1(1:num_eig));
V3=VV2(:, sort_idx2(2:num_eig+1));
delta=0.1; %cutoff value
if min(sorted_eigs1)<0
    id=find(abs(sorted_eigs1)>abs(min(sorted_eigs1)));
    sorted_eigs1=sorted_eigs1(id);
    V4=VV1(:, id);
    id=find(abs(sorted_eigs2)>abs(min(sorted_eigs2)));
    sorted_eigs2=sorted_eigs2(id);
    V3=VV2(:, id);
end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
t=ceil(log(delta)/log(min(abs(sorted_eigs1(1:num_eig)))))
t=0;
V4=repmat(sorted_eigs1(1:num_eig).', length(V4), 1).^t.*V4;
V3=repmat(sorted_eigs2(1:num_eig).', length(V3), 1).^t.*V3;
list=zeros((num_eig+1)*num_eig/2,2);
count=1;
for i=1:num_eig
    for j=i:num_eig
        list(count, 1)=i;
        list(count, 2)=j;
        if (j==i)
            c(count)=1;
        end;            
        count=count+1;
    end;
end;
r_VDM=conj(V4(:, list(:, 1))).*V4(:,list(:, 2));
r_DM=conj(V3(:, list(:, 1))).*V3(:, list(:, 2));
r_VDM(:, c==0)=r_VDM(:, c==0)*sqrt(2);
r_DM(:, c==0)=r_DM(:, c==0)*sqrt(2);
for i=1:M
    r_VDM(i, :)=r_VDM(i, :)/norm(r_VDM(i,:));
    r_DM(i, :)=r_DM(i, :)/norm(r_DM(i, :));
end;

end



