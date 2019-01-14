%Build three 3Nx3N matrices and do svd
function [rot,svals,s_out,svals2]= restoreRotations_eff_dev2(rr,cVec,N)
disp('Restoring rotations...');
tic

% curr_pool=gcp('nocreate');
% if isempty(curr_pool)
%     parpool('local',3);
% end
%% Fill 4D version of cMat for taking averages for Rii
npairs=nchoosek(N,2);
cMat5D=zeros(3,3,N,N,3);
cMat4D=zeros(3,3,npairs,3);
for i=1:N-1
    for j=i+1:N        
        ij=uppertri_ijtoind(i,j,N);
        cMat5D(:,:,i,j,cVec(3*ij-2))=rr(:,:,ij,1);
        cMat5D(:,:,i,j,cVec(3*ij-1))=rr(:,:,ij,2);
        cMat5D(:,:,i,j,cVec(3*ij))=rr(:,:,ij,3);
        cMat5D(:,:,j,i,cVec(3*ij-2))=rr(:,:,ij,1)';
        cMat5D(:,:,j,i,cVec(3*ij-1))=rr(:,:,ij,2)';
        cMat5D(:,:,j,i,cVec(3*ij))=rr(:,:,ij,3)';
        
        cMat4D(:,:,ij,cVec(3*ij-2))=rr(:,:,ij,1);
        cMat4D(:,:,ij,cVec(3*ij-1))=rr(:,:,ij,2);
        cMat4D(:,:,ij,cVec(3*ij))=rr(:,:,ij,3);
    end
end

%% Calculate Rii for color matrix
for i=1:N
    for c=1:3
        Rijs=squeeze(cMat5D(:,:,i,[1:i-1,i+1:N],c));
        Rii_est=multiprod(Rijs,permute(Rijs,[2,1,3]));
        Rii=sum(Rii_est,3)/(N-1);
        [U,~,~]=svd(Rii);
        cMat5D(:,:,i,i,c)=U(:,1)*U(:,1)';
    end
end
    
%% fill color matrix
cMat=zeros(3*N,3*N,3);
rot=zeros(3,3,N);
for i=1:N-1
    for j=i+1:N
        ij=uppertri_ijtoind(i,j,N);
        cMat(3*i-2:3*i,3*j-2:3*j,cVec(3*ij-2))=rr(:,:,ij,1);
        cMat(3*i-2:3*i,3*j-2:3*j,cVec(3*ij-1))=rr(:,:,ij,2);
        cMat(3*i-2:3*i,3*j-2:3*j,cVec(3*ij))=rr(:,:,ij,3);
    end
end
cMat(:,:,1)=cMat(:,:,1)+cMat(:,:,1)';
cMat(:,:,2)=cMat(:,:,2)+cMat(:,:,2)';
cMat(:,:,3)=cMat(:,:,3)+cMat(:,:,3)';

for c=1:3
    for i=1:N
        cMat(3*i-2:3*i,3*i-2:3*i,c)=cMat5D(:,:,i,i,c);
    end
end

%% Calc rank 1 matrices products of form Rik*Rkj for each row (color)
npairs=nchoosek(N,2);
signs=zeros(N,npairs,3);
%Rijs_prods_c=cell(N,N,c);
signs_c=cell(npairs,c);
pairs_idx_map=lin2sub_map(N);
parfor c=1:3
    pairs_idx_map_loc=pairs_idx_map;
    for p=1:npairs
        i=pairs_idx_map_loc(p,1);
        j=pairs_idx_map_loc(p,2);
        signs_c(p,c)={calcRijProds_ij(cMat5D,i,j,N,c)};
    end
end
%fill Rijs prods array after parfor      
for c=1:3
    for i=1:N-1
        for j=i+1:N
            ij=uppertri_ijtoind(i,j,N);
            signs([1:i-1,i+1:j-1,j+1:N],ij,c)=...
                cell2mat(signs_c(ij,c));
        end
    end
end

%Calculate relative signs between Rij and it's approximations Rik*Rkj
%s=size(Rijs_prods);
%signs=reshape(Rijs_prods,[9,s(3:5)]);
%signs=permute(sign(squeeze(sum(signs,1))),[2,1,3]);
%clear Rijs_prods
s=size(cMat4D);
est_signs=reshape(cMat4D,[9,s(3:4)]);
est_signs=sign(squeeze(sum(est_signs,1)));
signs=permute(signs,[2,1,3]);
%Sync signs to image 1 for testing
for c=1:3
    relative_signs=est_signs(:,c).*signs(:,:,c);
    signs(:,:,c)=relative_signs;
end

%sync_signs=zeros(N,N,N,3);
sync_block=repmat(1:N,N,1);
sync_signs2=repmat(sync_block,1,1,N,3);
clear sync_block
for c=1:3
    for r=1:N 
        %fill signs for sync with with each row/image #r
        i_idx=[1:r-1,r+1:N-1];
        for i=i_idx
            if i<=r
                j_idx=[i+1:r-1,r+1:N];
            else
                j_idx=i+1:N;
            end
            for j=j_idx
                ij=uppertri_ijtoind(i,j,N);
                %sync_signs(i,j,r,c)=signs(ij,r,c);
                %sync_signs(j,i,r,c)=signs(ij,r,c);
                sync_signs2(i,j,r,c)=j+0.5*(1-signs(ij,r,c))*N;
                sync_signs2(j,i,r,c)=i+0.5*(1-signs(ij,r,c))*N;
                %The function (1-x)/2 maps 1->0 and -1->1
            end
        end
    end
end
clear signs
%save('sync_signs2.mat','sync_signs2');
%sync_signs=sync_signs+repmat(eye(N),[1,1,N,3]); %DEBUG line
%% Choose signs for Rij according to sign estimates
% cMat5Dmp=gpuArray(single(cat(4,cMat5D,-cMat5D)));
% cMat_eff=gpuArray(zeros(3,3,N,N,'single'));
% rows_arr=gpuArray(zeros(3*N,N,3,'single'));

cMat5Dmp=cat(4,cMat5D,-cMat5D);
%cMat_eff=zeros(3,3,N,N);
rows_arr=zeros(3*N,N,3);
svals=zeros(N,2,3);
%cMat_arr=zeros(3*N,3*N,N,3);
%cMat=gpuArray(single(zeros(N,N,3)));

toc
curr_pool=gcp('nocreate');
if isempty(curr_pool)
    parpool('local',3);
end

parfor c=1:3 %color 
    for img=1:N
   
    %sync_signs_loc=sync_signs2;
     %image used for signs

%         for r=1:N %row of matrix to fill
%             cMat_eff(:,:,r,:)=cMat5Dmp(:,:,r,sync_signs2(r,:,img,c),c);
%         end
        cMat_eff=fill_signSync_matrix_c(cMat5Dmp,sync_signs2,img,c,N);
        %construct (3*N)x(3*N) rank 1 matrix 
        cMat_for_svd=zeros(3*N,3*N);%gpuArray(zeros(3*N,3*N,'single'));
%         for i=1:N-1
%             for j=i+1:N
%                 cMat_for_svd(3*i-2:3*i,3*j-2:3*j)=cMat_eff(:,:,i,j);
%             end
%         end
        for i=1:N
            row_3x3N=cMat_eff(:,:,i,:);
            row_3x3N=reshape(row_3x3N(:),3,3*N);
            cMat_for_svd(3*i-2:3*i,:)=row_3x3N;
        end
        cMat_for_svd=cMat_for_svd+cMat_for_svd';
        for i=1:N
            cMat_for_svd(3*i-2:3*i,3*i-2:3*i)=cMat_eff(:,:,i,i);
        end        
        [U,S,~]=svds(cMat_for_svd,2);
        svals(img,:,c)=gather(diag(S));
        rows_arr(:,img,c)=U(:,1); 
        
        %cMat_arr(:,:,img,c)=cMat_for_svd;%double(gather(cMat_for_svd));
    end
end
clear sync_signs2
%% Sync signs acoording to results for each image 
% delete(gcp('nocreate'));
% parpool('local',3);
%rows_arr=gather(rows_arr);
%sign_mat=zeros(npairs,npairs,3);
pairs_map=zeros(nchoosek(N,2),2*(N-2));
for i=1:N
    for j=i+1:N
        ij=uppertri_ijtoind(i,j,N);
        pairs_map(ij,:)=[uppertri_ijtoind(1:i-1,i,N),...
            uppertri_ijtoind(i,[i+1:j-1,j+1:N],N),...
            uppertri_ijtoind([1:i-1,i+1:j-1],j,N),...
            uppertri_ijtoind(j,j+1:N,N)];
    end
end

signs=zeros(npairs,3);
s_out=zeros(3,3);

% curr_pool=gcp('nocreate');
% if isempty(curr_pool)
%     parpool('local',3);
% end

parfor c=1:3
    sign_mat=zeros(npairs,2*(N-2));
    rows_arr_loc=rows_arr;
    pairs_map_loc=pairs_map;
    for i=1:N-1
        for j=i+1:N
            ij=uppertri_ijtoind(i,j,N);
            sij=rows_arr_loc(3*i-2:3*i,j,c);
            sji=rows_arr_loc(3*j-2:3*j,i,c);
            siks=rows_arr_loc(3*i-2:3*i,[1:i-1,i+1:j-1,j+1:N],c);
            sjks=rows_arr_loc(3*j-2:3*j,[1:i-1,i+1:j-1,j+1:N],c);
            sign_mat(ij,:)=[sign(sij'*siks),sign(sji'*sjks)];
        end
    end
    %sign_mat=sign_mat+sign_mat';%+eye(npairs);
    smat=@(x,flag) mult_smat_by_vec(x,sign_mat,pairs_map_loc,N);
    [U,S,~]=svds(smat,[npairs,npairs],3,'largest');
    signs(:,c)=U(:,1);
    s_out(:,c)=diag(S);
end

signs=sign(signs);

%% Assemble rotations 
for c=1:3
    idx=0;
    for i=1:N-1
        for j=i+1:N
            idx=idx+1;
            cMat(3*i-2:3*i,3*j-2:3*j,c)=...
                signs(idx,c)*cMat(3*i-2:3*i,3*j-2:3*j,c);
            cMat(3*j-2:3*j,3*i-2:3*i,c)=...
                signs(idx,c)*cMat(3*j-2:3*j,3*i-2:3*i,c);
        end
    end
end

%cMat_g=gpuArray(single(cMat));
[U1,S1,~]=svds(cMat(:,:,1),3);
[U2,S2,~]=svds(cMat(:,:,2),3);
[U3,S3,~]=svds(cMat(:,:,3),3);
svals2=zeros(3,3);
svals2(:,1)=diag(S1);
svals2(:,2)=diag(S2);
svals2(:,3)=diag(S3);
% U1=gather(U1);
% U2=gather(U2);
% U3=gather(U3);

for j=1:N
    rot(:,:,j)=[U1(3*j-2:3*j,1)'/norm(U1(3*j-2:3*j,1));...
        U2(3*j-2:3*j,1)'/norm(U2(3*j-2:3*j,1));...
        U3(3*j-2:3*j,1)'/norm(U3(3*j-2:3*j,1))];
    if det(rot(:,:,j))<0
        rot(3,:,j)=-rot(3,:,j);
    end
end
toc

%NOTE: Test if all metrices are rotations (det==1) 
%      If not, do SVD to find nearest orthogonal matrix

    return;
end

function [Rij]=calcRijProds_ij(cMat5D,i,j,N,c)

k_idx=[1:i-1,i+1:j-1,j+1:N];
Rik=squeeze(cMat5D(:,:,i,k_idx,c));
Rkj=squeeze(cMat5D(:,:,k_idx,j,c));
Rij=multiprod(Rik,Rkj);

%New code: In case we get a zero score
%Arbitrarily choose sign +1
ij_signs=squeeze(sum(reshape(Rij,9,N-2),1));
zeros_idx=(ij_signs==0);
n_zero_idx=sum(zeros_idx);
if sum(n_zero_idx)>0
    Rij(:,:,zeros_idx)=repmat(eye(3),n_zero_idx,1);
end

s=size(Rij);
Rij=reshape(Rij,9,s(3:end));
Rij=sign(sum(Rij,1));
end




