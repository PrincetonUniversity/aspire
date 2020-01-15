
function [out,Rijs_rows,color_perms,D,V]=syncColors(Rijs,K,nCPU)

disp('syncing colors...');
tic

% Run in parallel on nCPU workers. 
% curr_pool=gcp('nocreate');
% if isempty(curr_pool)
%     parpool('local',nCPU);
% end

%%  Generate array of one rank matrices from which we can extract rows. 
%   Matrices are of the form 0.5(Ri^TRj+Ri^TgkRj). Each such matrix can be
%   written in the form Qi^T*Ik*Qj where Ik is a 3x3 matrix with all zero 
%   entries except for the entry a_kk, k in {1,2,3}. 
N_pairs=nchoosek(K,2);
Rijs_rows=zeros(3,3,N_pairs,3);
for layer=1:3
    Rijs_rows(:,:,:,layer)=0.5*(Rijs(:,:,:,1)+Rijs(:,:,:,layer+1));
end

%%  Partition the set of matrices Rijs_rows into 3 sets of matrices, where 
%   each set there are only matrices Qi^T*Ik*Qj for a unique value of k in
%   {1,2,3}. 

%   First determine for each pair of tuples of the form {Qi^T*Ik*Qj} and 
%   {Qr^T*Il*Qj} where {i,j}\cap{r,l}==1, whether l==r. 
color_perms=match_colors_eff(Rijs_rows,K,nCPU);
%save('color_perms.mat','color_perms');
toc

%%  Compute eigenvectors of color matrix. This is just a matrix of dimensions
%   3(N over 2)x3(N over 2) where each entry corresponds to a pair of
%   matrices {Qi^T*Ir*Qj} and {Qr^T*Il*Qj} and eqauls \delta_rl. 
%   The 2 leading eigenvectors span a linear subspace which contains a
%   vector which encodes the partition. All the entries of the vector are
%   either 1,0 or -1, where the number encodes which the index r in Ir. 
%   This vector is a linear combination of the two leading eigen vectors,
%   and so we 'unmix' these vectors to retrieve it. 
disp('done matching colors , now calculating Eigen Vectors...');
neig=3;
cmat=@(x) mult_cmat_by_vec(color_perms,x,K);
[colors,D]=eigs(cmat,3*N_pairs,neig,'largestreal');
D=diag(D);
V=colors;
[cp,~]=unmixColors(colors(:,1:2),360);
toc
out=cp; %colors(:,whichVec);

    %  ***DEBUG CODE ***
%     n=size(color_perms,1);
%     color_perms_deb=zeros(n,3);
%     p_n=zeros(1,3);
%     for i=1:n
%         n=color_perms(i);
%         p_n(1)=floor(n/100);
%         p_n(3)=mod(n,10);
%         p_n(2)=(n-p_n(1)*100-p_n(3))/10;
%         color_perms_deb(i,:)=p_n;
%     end
%     
%     n=3*nchoosek(K,2);
%     color_mat=zeros(n,n);
%     block=-ones(3,3);
%     tperms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
%     cperms=reshape(color_perms_deb',3*nchoosek(K,3),1);
%     for i=1:K-2
%         for j=i+1:K-1
%             for k=j+1:K
%                 ij=uppertri_ijtoind(i,j,K);
%                 ik=uppertri_ijtoind(i,k,K);
%                 jk=uppertri_ijtoind(j,k,K);
%                 ijk=trip_idx(i,j,k,K);
%                 curr_block=block;
%                 curr_block(0+1,tperms(cperms(3*ijk-2),0+1))=1;
%                 curr_block(1+1,tperms(cperms(3*ijk-2),1+1))=1;
%                 curr_block(2+1,tperms(cperms(3*ijk-2),2+1))=1;
%                 color_mat(3*ij-2:3*ij,3*ik-2:3*ik)=curr_block;
%                 color_mat(3*ik-2:3*ik,3*ij-2:3*ij)=curr_block';
%                 
%                 curr_block=block;
%                 curr_block(0+1,tperms(cperms(3*ijk-1),0+1))=1;
%                 curr_block(1+1,tperms(cperms(3*ijk-1),1+1))=1;
%                 curr_block(2+1,tperms(cperms(3*ijk-1),2+1))=1;
%                 color_mat(3*ij-2:3*ij,3*jk-2:3*jk)=curr_block;
%                 color_mat(3*jk-2:3*jk,3*ij-2:3*ij)=curr_block';
%                 
%                 curr_block=block;
%                 curr_block(0+1,tperms(cperms(3*ijk-0),0+1))=1;
%                 curr_block(1+1,tperms(cperms(3*ijk-0),1+1))=1;
%                 curr_block(2+1,tperms(cperms(3*ijk-0),2+1))=1;
%                 color_mat(3*jk-2:3*jk,3*ik-2:3*ik)=curr_block;
%                 color_mat(3*ik-2:3*ik,3*jk-2:3*jk)=curr_block';
%                 
%             end
%         end
%     end
%     save('color_mat','color_mat');
%     [V,D]=eigs(color_mat,135);
%    D=diag(D);
%     ones_mat=ones(n,n);
%     [V2,D2]=eig(ones_mat-color_mat);
    %END DEBUG
%     V=V(:,1:2);
end

function [colors]=match_colors_eff(Rijs_rows,K,nworkers)

colors=zeros(nchoosek(K,3),1);
ntrip=nchoosek(K,3);
ntrip_per_worker=floor(ntrip/nworkers);
ntrip_last_worker=ntrip-ntrip_per_worker*(nworkers-1);
iter_lin_idx=zeros(nworkers+1,1);
for i=2:nworkers
    iter_lin_idx(i)=iter_lin_idx(i-1)+ntrip_per_worker;
end
iter_lin_idx(nworkers+1)=iter_lin_idx(nworkers)+ntrip_last_worker;

workers_res=cell(1,nworkers);

parfor i=1:nworkers
    lin_idx_loc=iter_lin_idx;
    workers_res{i}=...
        match_colors_c3_i(Rijs_rows,K,lin_idx_loc(i),lin_idx_loc(i+1));    
end

for i=1:nworkers
    colors(iter_lin_idx(i)+1:iter_lin_idx(i+1))=workers_res{i};
end

end

function [colors_i]=match_colors_c3_i(Rijs_rows,K,from,to)

Rijs_rows_tpd=permute(Rijs_rows,[2,1,3,4]);
trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
m=zeros(6,6);
colors_i=zeros(nchoosek(K,3),3);
votes=zeros(nchoosek(K,3),1);

ntrip=to-from;
trip_idx=lin2sub3_map(K);
trip_idx=trip_idx(from+1:to,:);

k1s=uppertri_ijtoind_vec(trip_idx(:,1),trip_idx(:,2),K);
k2s=uppertri_ijtoind_vec(trip_idx(:,1),trip_idx(:,3),K);
k3s=uppertri_ijtoind_vec(trip_idx(:,2),trip_idx(:,3),K);
ks=[k1s,k2s,k3s];

%   Compute relative color permutations. See Section 7.2 of paper. 
for t=1:ntrip
    
    k1=ks(t,1);
    k2=ks(t,2);
    k3=ks(t,3);
    
    %For r=1:3 compute 3*3 products v_{ji}(r)v_{ik}v_{kj}
    prod_arr=multiprod(permute(Rijs_rows(:,:,k2,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k3,:));
    prod_arr_tmp=prod_arr;
    prod_arr=multiprod(Rijs_rows_tpd(:,:,k1,:),reshape(prod_arr,3,3,9,1));
    prod_arr=permute(reshape(prod_arr,9,3,3,3),[1,4,2,3]);
    
    %Compare to v_{jj}(r)=v_{ji}v_{ij}
    self_prods=squeeze(multiprod(Rijs_rows_tpd(:,:,k1,:),Rijs_rows(:,:,k1,:)));
    self_prods=reshape(self_prods,9,3);
    prod_arr1=prod_arr;
    prod_arr1(:,1,:,:)=prod_arr1(:,1,:,:)-self_prods(:,1);
    prod_arr1(:,2,:,:)=prod_arr1(:,2,:,:)-self_prods(:,2);
    prod_arr1(:,3,:,:)=prod_arr1(:,3,:,:)-self_prods(:,3);
    norms1=squeeze(sum(prod_arr1.^2,1));
    prod_arr2=prod_arr;
    prod_arr2(:,1,:,:)=prod_arr2(:,1,:,:)+self_prods(:,1);
    prod_arr2(:,2,:,:)=prod_arr2(:,2,:,:)+self_prods(:,2);
    prod_arr2(:,3,:,:)=prod_arr2(:,3,:,:)+self_prods(:,3); 
    norms2=squeeze(sum(prod_arr2.^2,1));
    
    %Compare to v_{jj}(r)=v_{jk}v_{kj}
    self_prods=squeeze(multiprod(Rijs_rows(:,:,k3,:),Rijs_rows_tpd(:,:,k3,:)));
    self_prods=reshape(self_prods,9,3);
    prod_arr1=prod_arr;
    prod_arr1(:,:,:,1)=prod_arr1(:,:,:,1)-self_prods(:,1);
    prod_arr1(:,:,:,2)=prod_arr1(:,:,:,2)-self_prods(:,2);
    prod_arr1(:,:,:,3)=prod_arr1(:,:,:,3)-self_prods(:,3);
    norms1=norms1+squeeze(sum(prod_arr1.^2,1));
    prod_arr2=prod_arr;
    prod_arr2(:,:,:,1)=prod_arr2(:,:,:,1)+self_prods(:,1);
    prod_arr2(:,:,:,2)=prod_arr2(:,:,:,2)+self_prods(:,2);
    prod_arr2(:,:,:,3)=prod_arr2(:,:,:,3)+self_prods(:,3); 
    norms2=norms2+squeeze(sum(prod_arr2.^2,1));    
    
    %For r=1:3 compute 3*3 products v_{ij}(r)v_{jk}v_{ki} and compare to
    %Compare to v_{ii}(r)=v_{ij}v_{ji}
    prod_arr=permute(prod_arr_tmp,[2,1,3,4]);
    prod_arr=multiprod(Rijs_rows(:,:,k1,:),reshape(prod_arr,3,3,9,1));
    prod_arr=permute(reshape(prod_arr,9,3,3,3),[1,4,2,3]);
%     self_prods=squeeze(multiprod(Rijs_rows(:,:,k1,:),Rijs_rows_tpd(:,:,k1,:)));
%     self_prods=reshape(self_prods,9,3);
%     prod_arr1=prod_arr;
%     prod_arr1(:,1,:,:)=prod_arr1(:,1,:,:)-self_prods(:,1);
%     prod_arr1(:,2,:,:)=prod_arr1(:,2,:,:)-self_prods(:,2);
%     prod_arr1(:,3,:,:)=prod_arr1(:,3,:,:)-self_prods(:,3);
%     norms1=norms1+squeeze(sum(prod_arr1.^2,1));
%     prod_arr2=prod_arr;
%     prod_arr2(:,1,:,:)=prod_arr2(:,1,:,:)+self_prods(:,1);
%     prod_arr2(:,2,:,:)=prod_arr2(:,2,:,:)+self_prods(:,2);
%     prod_arr2(:,3,:,:)=prod_arr2(:,3,:,:)+self_prods(:,3); 
%     norms2=norms2+squeeze(sum(prod_arr2.^2,1));
    
    %Compare to v_{ii}(r)=v_{ik}v_{ki}
    self_prods=squeeze(multiprod(Rijs_rows(:,:,k2,:),Rijs_rows_tpd(:,:,k2,:)));
    self_prods=reshape(self_prods,9,3);
    prod_arr1=prod_arr;
    prod_arr1(:,:,1,:)=prod_arr1(:,:,1,:)-self_prods(:,1);
    prod_arr1(:,:,2,:)=prod_arr1(:,:,2,:)-self_prods(:,2);
    prod_arr1(:,:,3,:)=prod_arr1(:,:,3,:)-self_prods(:,3);
    norms1=norms1+squeeze(sum(prod_arr1.^2,1));
    prod_arr2=prod_arr;
    prod_arr2(:,:,1,:)=prod_arr2(:,:,1,:)+self_prods(:,1);
    prod_arr2(:,:,2,:)=prod_arr2(:,:,2,:)+self_prods(:,2);
    prod_arr2(:,:,3,:)=prod_arr2(:,:,3,:)+self_prods(:,3); 
    norms2=norms2+squeeze(sum(prod_arr2.^2,1));
    
    norms=min(norms1,norms2);
    
    for l=1:6
        p1=trip_perms(l,:);
        for r=1:6
            p2=trip_perms(r,:);
            m(l,r)=norms(1,p1(1),p2(1))+norms(2,p1(2),p2(2))+norms(3,p1(3),p2(3));
        end
    end
    idx=from+t; %trip_idx(i,j,k,K);
    [votes(idx),tmp]=min(m(:));
    %votes(idx)=sum(m(:));
    col=ceil(tmp/6);
    row=tmp-6*(col-1);
    colors_i(idx,1:2)=[100*row,10*col];
    
    % Calculate the relative permutation of Rik to Rij given 
    % by (sigma_ik)\circ(sigma_ij)^-1
    inv_jk_perm=inverse_perms(col,:);
    rel_perm=trip_perms(row,:);
    rel_perm=rel_perm(inv_jk_perm);
    colors_i(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
end
colors_i=colors_i(from+1:to,:);
colors_i=sum(colors_i,2);

end

%  Multiply color matrix by vector v "on the fly". 
%   Input:  v = vector to be multiplied by 'color matrix'. 
%           cperms = An (N over 3) vector. Each corresponds to a triplet of
%           indices i<j<k and indicates the relative permutation of tuples
%           {Qi^TIrQj},{Qj^TIlQk} and {Qk^TIlQi}. The color matrix can be 
%           completely reconstructed from
%           this information, which is used here to execute a single
%           multiplication of the matrix by the vector v instead of
%           explicitely computing and storing the prohibitively large
%           color matrix in memory. 
%           N = size of data. 
%    ***This loop intensive code should be written in C to gain a significant 
%       speed up.
function [out]=mult_cmat_by_vec(cperms,v,N)
    
    tperms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]-1; %  [1,2,3] permutations
    iperms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1]-1; %  Inverse permutations
    out=zeros(length(v),1);
    p_n=zeros(3,1);
    for i=1:N-2
        for j=i+1:N-1
            for k=j+1:N
                ijk=trip_idx(i,j,k,N);
                ij=3*uppertri_ijtoind_vec(i,j,N)-2;
                ik=3*uppertri_ijtoind_vec(i,k,N)-2;
                jk=3*uppertri_ijtoind_vec(j,k,N)-2;
                %   Extract permutation indexes from cperms
                n=cperms(ijk);
                p_n(1)=floor(n/100);
                p_n(3)=mod(n,10);
                p_n(2)=(n-p_n(1)*100-p_n(3))/10;
                
                %   Multiply vector by color matrix
                %   Upper triangular part
                p=tperms(p_n(1),:)+ik;
                out(ij)=out(ij)-v(p(2))-v(p(3))+v(p(1));
                out(ij+1)=out(ij+1)-v(p(1))-v(p(3))+v(p(2));
                out(ij+2)=out(ij+2)-v(p(1))-v(p(2))+v(p(3));
                p=tperms(p_n(2),:)+jk;
                out(ij)=out(ij)-v(p(2))-v(p(3))+v(p(1));
                out(ij+1)=out(ij+1)-v(p(1))-v(p(3))+v(p(2));
                out(ij+2)=out(ij+2)-v(p(1))-v(p(2))+v(p(3));
                p=iperms(p_n(3),:)+jk;
                out(ik)=out(ik)-v(p(2))-v(p(3))+v(p(1));
                out(ik+1)=out(ik+1)-v(p(1))-v(p(3))+v(p(2));
                out(ik+2)=out(ik+2)-v(p(1))-v(p(2))+v(p(3));
                %   Lower triangular part
                p=iperms(p_n(1),:)+ij;
                out(ik)=out(ik)-v(p(2))-v(p(3))+v(p(1));
                out(ik+1)=out(ik+1)-v(p(1))-v(p(3))+v(p(2));
                out(ik+2)=out(ik+2)-v(p(1))-v(p(2))+v(p(3));
                p=iperms(p_n(2),:)+ij;
                out(jk)=out(jk)-v(p(2))-v(p(3))+v(p(1));
                out(jk+1)=out(jk+1)-v(p(1))-v(p(3))+v(p(2));
                out(jk+2)=out(jk+2)-v(p(1))-v(p(2))+v(p(3));
                p=tperms(p_n(3),:)+ik;
                out(jk)=out(jk)-v(p(2))-v(p(3))+v(p(1));
                out(jk+1)=out(jk+1)-v(p(1))-v(p(3))+v(p(2));
                out(jk+2)=out(jk+2)-v(p(1))-v(p(2))+v(p(3));
                
            end
        end
    end
end

%%  The 'color vector' which partitions the rank 1 3x3 matrices into 3 sets
%   is one of 2 leading orthogonal eigenvectors of the color matrix. 
%   SVD retrieves two orthogonal linear combinations of these vectors which
%   can be 'unmixed' to retrieve the color vector by finding a suitable
%   2D rotation of these vectors (see Section 7.3 of D2 paper for details). 
function [colors,best_unmix]=unmixColors(ev,ntheta)

np=size(ev,1)/3;
dtheta=360/ntheta;
maxt=360/dtheta+1;
R_theta=@(theta)([cos(theta),-sin(theta);sin(theta),cos(theta)]);
s=inf;
scores=zeros(1,maxt);
idx=0;
for t=0:0.5:maxt
    idx=idx+1;
    
    unmix_ev=ev*R_theta(pi*t/180);
    s1=reshape(unmix_ev(:,1),3,np);
    [s1,p11]=sort(s1,1,'descend');
    s2=abs(reshape(unmix_ev(:,2),3,np));
    [s2,~]=sort(s2,1,'descend');
    score11=sum((s1(1,:)+s1(3,:)).^2+s1(2,:).^2);
    score12=sum((s2(1,:)-2*s2(2,:)).^2+(s2(1,:)-2*s2(3,:)).^2+...
        (s2(2,:)-s2(3,:)).^2); %Is this an error??? + instead of - in the first 2 members
    
    unmix_ev=ev*R_theta(pi*t/180);
    s1=abs(reshape(unmix_ev(:,1),3,np));
    [s1,~]=sort(s1,1,'descend');
    s2=reshape(unmix_ev(:,2),3,np);
    [s2,p22]=sort(s2,1,'descend');    
    score21=sum((s2(1,:)+s2(3,:)).^2+s2(2,:).^2);
    score22=sum((s1(1,:)-2*s1(2,:)).^2+(s1(1,:)-2*s1(3,:)).^2+...
        (s1(2,:)-s1(3,:)).^2);
    
    [scores(idx),whichVec]=min([score11+score12,score21+score22]);
    if scores(idx)<s
        s=scores(idx);
        if whichVec==1
            p=p11;
        else
            p=p22;
        end
        best_unmix=unmix_ev(:,whichVec);
    end    
end

%Assign integers between 1:3 to permutations
colors=zeros(3,np);
for i=1:np
    p_i=p(:,i);
    p_i_sqr=p_i(p_i);
    if sum((p_i_sqr-[1;2;3]).^2)==0 % non cyclic pemutation
        colors(:,i)=p_i;
    else
        colors(:,i)=p_i_sqr;
    end
end
colors=colors(:);
end