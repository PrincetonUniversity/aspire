function  [ r_VDM, r_DM, V4 ]=VDM_LP( H1, num_eig)
%VDM and DM embedding
%   Input:
%       H1: weight matrix with orientation H1_ij = H_ij e^{i\alpha_ij}.
%       num_eig: the number of eigenvectors to use.
%   Output:
%       r_VDM: the embedding from VDM
%       r_DM: the embedding from DM
%       V4: the top eigenvectors of size N x num_eig.
%   Zhizhen Zhao June 01 2013

%% Eigenvector calculation
H=abs(H1);
M=size(H, 1);
C=sum(H, 2);
D=sparse(1:M, 1:M, 1./sqrt(C), M, M);
H=D*H*D;
H1=D*H1*D;
[VV1, D1]=eigs(H1, num_eig);
[VV2, D2]=eigs(H, num_eig+1);
eig_vals1=diag(D1);
eig_vals2=diag(D2);
[sorted_eigs1, sort_idx1] = sort(real(eig_vals1),'descend');
[sorted_eigs2, sort_idx2] = sort(real(eig_vals2),'descend');
V4=VV1(:, sort_idx1(1:num_eig));
V3=VV2(:, sort_idx2(2:num_eig+1));

if min(sorted_eigs1)<0
    id=find(abs(sorted_eigs1)>abs(min(sorted_eigs1)));
    sorted_eigs1=sorted_eigs1(id);
    V4=VV1(:, id);
    id=find(abs(sorted_eigs2)>abs(min(sorted_eigs2)));
    sorted_eigs2=sorted_eigs2(id);
    V3=VV2(:, id);
end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
delta=0.1; %cutoff value
t=ceil(log(delta)/log(min(abs(sorted_eigs1(1:num_eig)))));
V4=repmat(sorted_eigs1(1:num_eig).', length(V4), 1).^t.*V4;
V3=repmat(sorted_eigs2(2:num_eig+1).', length(V3), 1).^t.*V3;
list=zeros((num_eig+1)*num_eig/2,2);
c=zeros((num_eig+1)*num_eig/2, 1);
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
r_DM=V3;
r_VDM(:, c==0)=r_VDM(:, c==0)*sqrt(2);
for i=1:M
    r_VDM(i, :)=r_VDM(i, :)/norm(r_VDM(i,:));
    r_DM(i, :)=r_DM(i, :)/norm(r_DM(i, :));
end;

%normlize rows of the eigenvectors (weighted by eigenvalues).
% for i=1:M
%     V4(i, :)=V4(i, :)/norm(V4(i, :));
% end;

end



