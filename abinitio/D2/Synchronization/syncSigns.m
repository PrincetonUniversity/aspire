%%  This function executes the final stage of the algorithm, Signs 
%   synchroniztion. At the end all rows of the rotations Ri are exctracted
%   and the matrices Ri are assembled. 
   
function [rot,svals,s_out,svals2]= syncSigns(rr,cVec,N)
log_message('Synchronizing signs...');
tic
%%  Partition the union of tuples {0.5*(Ri^TRj+Ri^TgkRj), k=1:3} according 
%   to the color partition established in color synchroniztion procedure. 
%   The partition is stored in two different arrays each with the purpose
%   of a computational speed up for two different computations performed 
%   later (space considerations are of little concern since arrays are ~ 
%   o(N^2) which doesn't pose a constraint for inputs on the scale of 10^3-10^4. 
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

%%  Compute estimates for the tuples {0.5*(Ri^TRi+Ri^TgkRi), k=1:3} for 
%   i=1:N. For 1<=i,j<=N and c=1,2,3 write Qij^c=0.5*(Ri^TRj+Ri^TgmRj).  
%   For each i in {1:N} and each k in {1,2,3} the estimator is the
%   average over all j~=i of Qij^c*(Qij^c)^T. 
%   Since in practice the result of the average is not really rank 1, we
%   compute the best rank approximation to this average. 
for i=1:N
    for c=1:3
        Rijs=squeeze(cMat5D(:,:,i,[1:i-1,i+1:N],c));
        Rii_est=multiprod(Rijs,permute(Rijs,[2,1,3]));
        Rii=sum(Rii_est,3)/(N-1);
        [U,~,~]=svd(Rii);
        cMat5D(:,:,i,i,c)=U(:,1)*U(:,1)';
    end
end
    
%%  Construct the 3Nx3N row synchroniztion matrices (as done for C_2), one
%   for all first rows of the matrices Ri, one for all second rows and one 
%   for all third rows. The ij'th block of the k'th matrix is Qij^c. 
%   In C_2 one such matrix is constructed for the 3rd rows
%   and is rank 1 by construction. In practice, thus far, for each c and
%   (i,j) we either have Qij^c or -Qij^c independently. 
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

%%  To decompose cMat as a rank 1 matrix we need to adjust the signs of the 
%   Qij^c so that sign(Qij^c*Qjk^c) = sign(Qik^c) for all c=1,2,3 and (i,j).
%   In practice we compare the sign of the sum of the entries of Qij^c*Qjk^c 
%   to the sum of entries of Qik^c. 
npairs=nchoosek(N,2);
signs=zeros(N,npairs,3);
signs_c=cell(npairs,c);
pairs_idx_map=lin2sub_map(N);
%   First we compute the signs for Qij^c*Qjk^c. 
parfor c=1:3
    pairs_idx_map_loc=pairs_idx_map;
    for p=1:npairs
        i=pairs_idx_map_loc(p,1);
        j=pairs_idx_map_loc(p,2);
        signs_c(p,c)={calcRijProds_ij(cMat5D,i,j,N,c)};
    end
end
%	For computational comfort the signs for each c=1,2,3 are stored in a
%	Nx(N over 2) array, where the ij'th column corresponds to the signs of
%   Qij^c * Qjk^c for k~=i,j. The entries in the k=i,j rows of the ij'th 
%   column are zero, the value zero is arbitrary, since these entries are
%   not used by the algorithm, and only exist for comfort (of storage and
%   excess). 
for c=1:3
    for i=1:N-1
        for j=i+1:N
            ij=uppertri_ijtoind(i,j,N);
            signs([1:i-1,i+1:j-1,j+1:N],ij,c)=...
                cell2mat(signs_c(ij,c));
        end
    end
end

%   Now compute the signs of Qij^c.
s=size(cMat4D);
est_signs=reshape(cMat4D,[9,s(3:4)]);
est_signs=sign(squeeze(sum(est_signs,1)));
signs=permute(signs,[2,1,3]);

%%	Compute relative signs between Qik^c and Qij^c*Qjk^c. 
log_message('Computing relative signs...');
%   c=1,2,3. 
for c=1:3
%     relative_signs=est_signs(:,c).*signs(:,:,c);
%     signs(:,:,c)=relative_signs;
    signs(:,:,c) = est_signs(:,c).*signs(:,:,c);
end

%   Qik^c can be compared with Qir^c*Qrk^c for each r~=i,k, that is,
%   N-2 options. Another way to look at this, is that the r'th image
%   participates in all comparisons of the form sign(Qir^c*Qrk^c)~sign(Qik)
%   for r~=i,k for each c=1,2,3 (see Section 8 in D2 paper). 
%   For each image r construct a 3Nx3N matrix. If
%   sign(Qir^c*Qrk^c)~sign(Qik)=1, its ik'th 3x3 block is set to Qik,
%   otherwise, it is set to -Qik. 
sync_block=repmat(1:N,N,1);
sync_signs2=repmat(sync_block,1,1,N,3);
clear sync_block
for c=1:3
    for r=1:N 
        %   Fill signs for synchroniztion for the r'th image.
        %   Go over all i,j~=r. 
        i_idx=[1:r-1,r+1:N-1];   
        for i=i_idx %   i~=r
            if i<=r %   j~=r,i
                j_idx=[i+1:r-1,r+1:N];
            else
                j_idx=i+1:N;
            end
            for j=j_idx
                ij=uppertri_ijtoind(i,j,N);
                sync_signs2(i,j,r,c)=j+0.5*(1-signs(ij,r,c))*N;
                sync_signs2(j,i,r,c)=i+0.5*(1-signs(ij,r,c))*N;
                %The function (1-x)/2 maps 1->0 and -1->1
            end
        end
    end
end
clear signs

cMat5Dmp=cat(4,cMat5D,-cMat5D);
rows_arr=zeros(3*N,N,3);
svals=zeros(N,2,3);

% curr_pool=gcp('nocreate');
% if isempty(curr_pool)
%     parpool('local',3);
% end

log_message('Constructing and decomposing N sign synchroniztion matrices...');
parfor c=1:3 %color 
    for r=1:N
   
        %   Image r used for signs
        cMat_eff=fill_signSync_matrix_c(cMat5Dmp,sync_signs2,r,c,N);
        %   Construct (3*N)x(3*N) rank 1 matrices from Qik
        cMat_for_svd=zeros(3*N,3*N);%gpuArray(zeros(3*N,3*N,'single'));

        for i=1:N
            row_3x3N=cMat_eff(:,:,i,:);
            row_3x3N=reshape(row_3x3N(:),3,3*N);
            cMat_for_svd(3*i-2:3*i,:)=row_3x3N;
        end
        cMat_for_svd=cMat_for_svd+cMat_for_svd';
        %   Extract leading eigenvector of rank 1 matrix. For each r and c 
        %   this gives an estimate for the c'th row of the rotation Rr, up
        %   to sign +/-.
        for i=1:N
            cMat_for_svd(3*i-2:3*i,3*i-2:3*i)=cMat_eff(:,:,i,i);
        end        
        [U,S,~]=svds(cMat_for_svd,2);
        svals(r,:,c)=gather(diag(S));
        rows_arr(:,r,c)=U(:,1); 
    end
end
clear sync_signs2
%%  Sync signs acoording to results for each image. Dot products between 
%   signed row estimates are used to construct an (N over 2)x(N over 2) 
%   sign synchronization matrix S. If (v_i)k and (v_j)k are the i'th and 
%   j'th estimates for the c'th row of Rk, then the entry (i,k),(k,j) entry 
%   of S is <(v_i)k,(v_j)k>, where the rows and columns of S are indexed by
%   double indexes (i,j), 1<=i<j<=(N over 2). 
pairs_map=zeros(nchoosek(N,2),2*(N-2));
for i=1:N
    for j=i+1:N
        ij=uppertri_ijtoind(i,j,N);
        pairs_map(ij,:)=[uppertri_ijtoind_vec(1:i-1,i,N),...
            uppertri_ijtoind_vec(i,[i+1:j-1,j+1:N],N),...
            uppertri_ijtoind_vec([1:i-1,i+1:j-1],j,N),...
            uppertri_ijtoind_vec(j,j+1:N,N)];
    end
end

signs=zeros(npairs,3);
s_out=zeros(3,3);

log_message('Constructing and decomposing 3 sign synchroniztion matrices...');
%   The matrix S requires space on order of O(N^4). Instead of storing it
%   in memory we compute its SVD using the function smat which multiplies
%   (N over 2)x1 vectors by S.
parfor c=1:3
    %   Preapare data for smat to act on vectors. 
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
    %   smat acts using the auxiliary function mult_smat_by_vec. 
    smat=@(x,flag) mult_smat_by_vec(x,sign_mat,pairs_map_loc,N);
    [U,S,~]=svds(smat,[npairs,npairs],3,'largest');
    signs(:,c)=U(:,1);
    s_out(:,c)=diag(S);
end
signs=sign(signs);

%% Adjust the signs of Qij^c in the matrices cMat(:,:,c) for all c=1,2,3 
%  and 1<=i<j<=N according to the results of the signs from the last stage.  
log_message('Constructing and decomposing 3 row synchroniztion matrices...');
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

%%  cMat(:,:,c) are now rank 1. Decompose using SVD and take leading eigenvector. 
[U1,S1,~]=svds(cMat(:,:,1),3);
[U2,S2,~]=svds(cMat(:,:,2),3);
[U3,S3,~]=svds(cMat(:,:,3),3);
svals2=zeros(3,3);
svals2(:,1)=diag(S1);
svals2(:,2)=diag(S2);
svals2(:,3)=diag(S3);

%%  The c'th row of the rotation Rj is Uc(3*j-2:3*j,1)/norm(Uc(3*j-2:3*j,1)),
%   (Rows must be normalized to length 1). 
log_message('Assembeling rows to rotations matrices...');
for j=1:N
    rot(:,:,j)=[U1(3*j-2:3*j,1)'/norm(U1(3*j-2:3*j,1));...
        U2(3*j-2:3*j,1)'/norm(U2(3*j-2:3*j,1));...
        U3(3*j-2:3*j,1)'/norm(U3(3*j-2:3*j,1))];
    if det(rot(:,:,j))<0
        rot(3,:,j)=-rot(3,:,j);
    end
end
toc
end
%%  Auxiliary functions
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

function [cMat_eff]=fill_signSync_matrix_c(cMat5Dmp,sync_signs2,img,c,N)
cMat_eff=zeros(3,3,N,N);
for r=1:N %row of matrix to fill
    cMat_eff(:,:,r,:)=cMat5Dmp(:,:,r,sync_signs2(r,:,img,c),c);
end
end

% mult_smat_by_vec multiplies the signs sync matrix to a vector using: 
% Input: rows_arr= an array of size 3NxN where rows_arr(3*i-2:3*i,j) is the
%        j'th estimate for the c'th row of the matrix Ri. 
%        pairs_map= array of size NxN where
%        pairs_map(i,j)==uppertri_ijtoind(i,j)

function [v_out]=mult_smat_by_vec(v,sign_mat,pairs_map,N)

v_out=zeros(size(v));

for i=1:N
    for j=i+1:N    
       ij=uppertri_ijtoind(i,j,N);
       v_out(ij)=sign_mat(ij,:)*v(pairs_map(ij,:));
    end
end
end


