
function [out,Rijs_rows,color_perms,D,V,unmix_colors_all]=sync_colors_unmix_new2(Rijs,K,scheme)
    disp('syncing colors...');
    tic
    %generate array of one rank matrices from which we can extract rows
    N_pairs=nchoosek(K,2);
    Rijs_rows=zeros(3,3,N_pairs,3);
    for layer=1:3
        Rijs_rows(:,:,:,layer)=0.5*(Rijs(:,:,:,1)+Rijs(:,:,:,layer+1));
    end
    
    curr_pool=gcp('nocreate');
    if isempty(curr_pool)
        parpool('local',12);
    end
    
    %match colors for 3 color scheme
    switch scheme
        case 1
            color_perms=match_colors_c1(Rijs_rows,K);
        case 2
            color_perms=match_colors_c2(Rijs_rows,K);
        case 3
            %color_perms=match_colors_c3(Rijs_rows,K);
            color_perms=match_colors_eff(Rijs_rows,K,12);
            save('color_perms.mat','color_perms');
            %color_perms=parallel_match_colors_c3(Rijs_rows,K);
            disp('done matching colors');
            toc
            disp('Calculating Eigen Vectors...');
            tic     
            
            %DEBUG
            %color_perms=repmat(111,120,1);
            %K=10;            
            neig=3;
            cmat=@(x) mult_cmat_by_vec_tmp(color_perms,x,K);
            [colors,D]=eigs(cmat,3*N_pairs,neig,'largestreal');
            D=diag(D);
            V=colors;
            [cp,~]=unmixColors_new1(colors(:,1:2),360);
            toc
            out=cp; %colors(:,whichVec);
            return;
        case 4
            load('colors_500_snr_0.5.mat');
            %[colors1,angles]=recolor(colors);
            var3=0; var4=0;       
            for i=1:N_pairs
                [sorted,classes]=sort_3(colors(3*i-2:3*i,1));
                colors(3*i-2:3*i,1)=classes;
                var3=var3+sorted(2);%abs(sorted(2)-0.5*(sorted(1)+sorted(3)));
                [sorted,classes]=sort_3(colors(3*i-2:3*i,2));
                colors(3*i-2:3*i,2)=classes;
                var4=var4+sorted(2); %abs(sorted(2)-0.5*(sorted(1)+sorted(3)));
            end
            if var3<var4
                out=colors(:,1);
            else
                out=colors(:,2);
            end
            color_perms=[];
            return;
    end
    %create a sparse matrix for color mat. There are nchoosek(K,3) triplets
    %i,j,k for each there are 2*3 blocks of 3x3 with 6 non-zero entiries
%     sparse_ipt=ones(6*2*3*nchoosek(K,3),3);
%     p=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
%     cp=color_perms;
%     
%     for i=1:K-2
%         for j=i+1:K-1
%             for k=j+1:K
%                 ij=uppertri_ijtoind(i,j,K);
%                 ik=uppertri_ijtoind(i,k,K);
%                 jk=uppertri_ijtoind(j,k,K);
%                 ijk=trip_idx(i,j,k,K);
%                 p_ij_ik=p(cp(ijk,1));
%                 sparse_ipt(18*idx-17,:)=[3*ij-2,3*ik-3+p_ij_ik(1),0];
%                 sparse_ipt(18*idx-16,:)=[3*ij-1,3*ik-3+p_ij_ik(2),0];
%                 sparse_ipt(18*idx-15,:)=[3*ij,3*ik-3+p_ij_ik(3),0];
%                 sparse_ipt(18*idx-14,:)=[3*ik-3+p_ij_ik(1),3*ij-2,0];
%                 sparse_ipt(18*idx-13,:)=[3*ik-3+p_ij_ik(2),3*ij-1,0];
%                 sparse_ipt(18*idx-12,:)=[3*ik-3+p_ij_ik(3),3*ij,0];
%                 
%                 p_ij_jk=p(cp(ijk,2));
%                 sparse_ipt(18*idx-11,:)=[3*ij-2,3*jk-3+p_ij_jk(1),0];
%                 sparse_ipt(18*idx-10,:)=[3*ij-1,3*jk-3+p_ij_jk(2),0];
%                 sparse_ipt(18*idx-9,:)=[3*ij,3*jk-3+p_ij_jk(3),0];
%                 sparse_ipt(18*idx-8,:)=[3*jk-3+p_ij_jk(1),3*ij-2,0];
%                 sparse_ipt(18*idx-7,:)=[3*jk-3+p_ij_jk(2),3*ij-1,0];
%                 sparse_ipt(18*idx-6,:)=[3*jk-3+p_ij_jk(3),3*ij,0];
%                 
%                 p_jk_ik=p(cp(ijk,3));
%                 sparse_ipt(18*idx-5,:)=[3*jk-2,3*ik-3+p_jk_ik(1),0];
%                 sparse_ipt(18*idx-4,:)=[3*jk-1,3*ik-3+p_jk_ik(2),0];
%                 sparse_ipt(18*idx-3,:)=[3*jk,3*ik-3+p_jk_ik(3),0];
%                 sparse_ipt(18*idx-2,:)=[3*ik-3+p_jk_ik(1),3*jk-2,0];
%                 sparse_ipt(18*idx-1,:)=[3*ik-3+p_jk_ik(2),3*jk-1,0];
%                 sparse_ipt(18*idx-0,:)=[3*ik-3+p_jk_ik(3),3*jk,0];
%             end
%         end
%     end
%     S=sparse(sparse_ipt(:,1),sparse_ipt(:,2),sparse_ipt(:,3));
%     opts=struct('issym',1,'isreal',1);
%     [V,D,flag]=eigs(S,3,'lm',opts);
        
%     in=reshape(color_perms',3*nchoosek(K,3),1);
%     colors =colors_power_method (K,in-1,3,0);
%     
%     N_pairs=nchoosek(K,2);
%     var1=0; var2=0; 
%     for i=1:N_pairs
%         [sorted,classes]=sort_3(colors(3*i-2:3*i,1));
%         colors(3*i-2:3*i,1)=classes;
%         var1=var1+abs(sorted(2)-0.5*(sorted(1)+sorted(3)));
%         [sorted,classes]=sort_3(colors(3*i-2:3*i,2));
%         colors(3*i-2:3*i,2)=classes;
%         var2=var2+abs(sorted(2)-0.5*(sorted(1)+sorted(3)));
%     end
%     
%     if var1<var2
%         colors=colors(:,1);
%     else
%         colors=colors(:,2);
%     end
    
    %DEBUG CODE
    n=size(color_perms,1);
    color_perms_deb=zeros(n,3);
    p_n=zeros(1,3);
    for i=1:n
        n=color_perms(i);
        p_n(1)=floor(n/100);
        p_n(3)=mod(n,10);
        p_n(2)=(n-p_n(1)*100-p_n(3))/10;
        color_perms_deb(i,:)=p_n;
    end
    
    n=3*nchoosek(K,2);
    color_mat=zeros(n,n);
    block=-ones(3,3);
    tperms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
    cperms=reshape(color_perms_deb',3*nchoosek(K,3),1);
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                ij=uppertri_ijtoind(i,j,K);
                ik=uppertri_ijtoind(i,k,K);
                jk=uppertri_ijtoind(j,k,K);
                ijk=trip_idx(i,j,k,K);
                curr_block=block;
                curr_block(0+1,tperms(cperms(3*ijk-2),0+1))=1;
                curr_block(1+1,tperms(cperms(3*ijk-2),1+1))=1;
                curr_block(2+1,tperms(cperms(3*ijk-2),2+1))=1;
                color_mat(3*ij-2:3*ij,3*ik-2:3*ik)=curr_block;
                color_mat(3*ik-2:3*ik,3*ij-2:3*ij)=curr_block';
                
                curr_block=block;
                curr_block(0+1,tperms(cperms(3*ijk-1),0+1))=1;
                curr_block(1+1,tperms(cperms(3*ijk-1),1+1))=1;
                curr_block(2+1,tperms(cperms(3*ijk-1),2+1))=1;
                color_mat(3*ij-2:3*ij,3*jk-2:3*jk)=curr_block;
                color_mat(3*jk-2:3*jk,3*ij-2:3*ij)=curr_block';
                
                curr_block=block;
                curr_block(0+1,tperms(cperms(3*ijk-0),0+1))=1;
                curr_block(1+1,tperms(cperms(3*ijk-0),1+1))=1;
                curr_block(2+1,tperms(cperms(3*ijk-0),2+1))=1;
                color_mat(3*jk-2:3*jk,3*ik-2:3*ik)=curr_block;
                color_mat(3*ik-2:3*ik,3*jk-2:3*jk)=curr_block';
                
            end
        end
    end
    save('color_mat','color_mat');
    [V,D]=eigs(color_mat,135);
%    D=diag(D);
%     ones_mat=ones(n,n);
%     [V2,D2]=eig(ones_mat-color_mat);
    %END DEBUG
    V=V(:,1:2);
end

%criterion: product of all 3 colors is maximal
function [colors]=match_colors_c1(Rijs_rows,K)
    tic
    Rijs_rows_tpd=multi_transpose(Rijs_rows);
    trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
    inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
    m=zeros(6,6);
    colors=zeros(nchoosek(K,3),3);
    votes=zeros(nchoosek(K,3),1);
    
    %we test all 6x6 permutations of Rik{1,2,3) and Rjk{1,2,3} vs.
    %Rij{1,2,3} and choose the best one (maximal sum on norms of triple
    %products of form sum(Rij_a*Rik_b*Rjk_c) where a+b+c=6, a,b,c are ints
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
%                 prod_arr=multiprod(permute(Rijs_rows(:,:,k3,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k2,:));
%                 prod_arr=multiprod(Rijs_rows(:,:,k1,:),reshape(prod_arr,3,3,9,1));
%                 prod_arr=permute(reshape(prod_arr,9,3,3,3),[1 4 2 3]);
                A=permute(Rijs_rows(:,:,k2,:),[1,2,4,3]);
                B=Rijs_rows_tpd(:,:,k3,:);
                prod_arr=multiprod(permute(Rijs_rows(:,:,k2,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k3,:));
                C=prod_arr;
                C_vec=reshape(prod_arr,3,3,9,1);
                prod_arr=multiprod(Rijs_rows_tpd(:,:,k1,:),reshape(prod_arr,3,3,9,1));
                D=prod_arr;
                prod_arr=permute(reshape(prod_arr,9,3,3,3),[1,4,2,3]);
                
                norms=squeeze(sum(prod_arr.^2,1));
                for l=1:6
                    p1=trip_perms(l,:);
                    for r=1:6
                        p2=trip_perms(r,:);
                        m(l,r)=norms(1,p1(1),p2(1))+norms(2,p1(2),p2(2))+norms(3,p1(3),p2(3));
                    end
                end 
                idx=trip_idx(i,j,k,K);
                [votes(idx),tmp]=max(m(:));
                votes(idx)=sum(m(:));
                col=ceil(tmp/6);
                row=tmp-6*(col-1);
                colors(idx,1:2)=[row,col];
                
                %Need to calculate the relative permutation of Rik to Rij
                %by (sigma_ik)\circ(sigma_ij)^-1 
                inv_jk_perm=inverse_perms(col,:);
                rel_perm=trip_perms(row,:);
                rel_perm=rel_perm(inv_jk_perm);
                colors(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
            end
        end
    end
    toc
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

Rijs_rows_tpd=permute(Rijs_rows,[2,1,3,4]);%multi_transpose(Rijs_rows);
trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
m=zeros(6,6);
colors_i=zeros(nchoosek(K,3),3);
votes=zeros(nchoosek(K,3),1);

ntrip=to-from;
trip_idx=lin2sub3_map(K);
trip_idx=trip_idx(from+1:to,:);

k1s=uppertri_ijtoind(trip_idx(:,1),trip_idx(:,2),K);
k2s=uppertri_ijtoind(trip_idx(:,1),trip_idx(:,3),K);
k3s=uppertri_ijtoind(trip_idx(:,2),trip_idx(:,3),K);
ks=[k1s,k2s,k3s];

%we test all 6x6 permutations of Rik{1,2,3) and Rjk{1,2,3} vs.
%Rij{1,2,3} and choose the best one (maximal sum on norms of triple
%products of form sum(Rij_a*Rik_b*Rjk_c) where a+b+c=6, a,b,c are ints
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
    
    %Need to calculate the relative permutation of Rik to Rij
    %by (sigma_ik)\circ(sigma_ij)^-1
    inv_jk_perm=inverse_perms(col,:);
    rel_perm=trip_perms(row,:);
    rel_perm=rel_perm(inv_jk_perm);
    colors_i(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
end
colors_i=colors_i(from+1:to,:);
colors_i=sum(colors_i,2);

end

function [colors]=match_colors_c3(Rijs_rows,K)

    Rijs_rows_tpd=multi_transpose(Rijs_rows);
    trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
    inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
    m=zeros(6,6);
    colors=zeros(nchoosek(K,3),3);
    votes=zeros(nchoosek(K,3),1);
    
    %we test all 6x6 permutations of Rik{1,2,3) and Rjk{1,2,3} vs.
    %Rij{1,2,3} and choose the best one (maximal sum on norms of triple
    %products of form sum(Rij_a*Rik_b*Rjk_c) where a+b+c=6, a,b,c are ints
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
%                 prod_arr=multiprod(permute(Rijs_rows(:,:,k3,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k2,:));
%                 prod_arr=multiprod(Rijs_rows(:,:,k1,:),reshape(prod_arr,3,3,9,1));
%                 prod_arr=permute(reshape(prod_arr,9,3,3,3),[1 4 2 3]);
                A=permute(Rijs_rows(:,:,k2,:),[1,2,4,3]);
                B=Rijs_rows_tpd(:,:,k3,:);
                prod_arr=multiprod(permute(Rijs_rows(:,:,k2,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k3,:));
                C=prod_arr;
                C_vec=reshape(prod_arr,3,3,9,1);
                prod_arr=multiprod(Rijs_rows_tpd(:,:,k1,:),reshape(prod_arr,3,3,9,1));
                D=prod_arr;
                prod_arr=permute(reshape(prod_arr,9,3,3,3),[1,4,2,3]);
                
                norms=squeeze(sum(prod_arr.^2,1));
                for l=1:6
                    p1=trip_perms(l,:);
                    for r=1:6
                        p2=trip_perms(r,:);
                        m(l,r)=norms(1,p1(1),p2(1))+norms(2,p1(2),p2(2))+norms(3,p1(3),p2(3));
                    end
                end 
                idx=trip_idx(i,j,k,K);
                [votes(idx),tmp]=max(m(:));
                votes(idx)=sum(m(:));
                col=ceil(tmp/6);
                row=tmp-6*(col-1);
                colors(idx,1:2)=[100*row,10*col];
                
                %Need to calculate the relative permutation of Rik to Rij
                %by (sigma_ik)\circ(sigma_ij)^-1 
                inv_jk_perm=inverse_perms(col,:);
                rel_perm=trip_perms(row,:);
                rel_perm=rel_perm(inv_jk_perm);
                colors(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
            end
        end
    end
    colors=sum(colors,2);

end

function [colors]=parallel_match_colors_c3(Rijs_rows,K)

    
%     colors=zeros(nchoosek(K,3),3);
%     votes=zeros(nchoosek(K,3),1);
    colors_c=cell(K,K,K);
    
    %we test all 6x6 permutations of Rik{1,2,3) and Rjk{1,2,3} vs.
    %Rij{1,2,3} and choose the best one (maximal sum on norms of triple
    %products of form sum(Rij_a*Rik_b*Rjk_c) where a+b+c=6, a,b,c are ints
    parfor i=1:K-2
        Rijs_rows_tpd=multi_transpose(Rijs_rows);
        trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
        inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
        m=zeros(6,6);
        for j=2:K-1
            if j<=i
                continue;
            end
            for k=3:K
                if k<=j
                    continue;
                end
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                prod_arr=multiprod(permute(Rijs_rows(:,:,k2,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k3,:));
                prod_arr=multiprod(Rijs_rows_tpd(:,:,k1,:),reshape(prod_arr,3,3,9,1));
                prod_arr=permute(reshape(prod_arr,9,3,3,3),[1,4,2,3]);
                
                norms=squeeze(sum(prod_arr.^2,1));
                for l=1:6
                    p1=trip_perms(l,:);
                    for r=1:6
                        p2=trip_perms(r,:);
                        m(l,r)=norms(1,p1(1),p2(1))+norms(2,p1(2),p2(2))+norms(3,p1(3),p2(3));
                    end
                end 
                %idx=trip_idx(i,j,k,K);
                [~,tmp]=max(m(:));
%                 votes(idx)=sum(m(:));
                col=ceil(tmp/6);
                row=tmp-6*(col-1);
                colors_out=[100*row,10*col,0];
                %colors(idx,1:2)=[100*row,10*col];
                
                %Need to calculate the relative permutation of Rik to Rij
                %by (sigma_ik)\circ(sigma_ij)^-1 
                inv_jk_perm=inverse_perms(col,:);
                rel_perm=trip_perms(row,:);
                rel_perm=rel_perm(inv_jk_perm);
                colors_out(3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
                %colors(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
                colors_c(i,j,k)={colors_out};
            end
        end
    end
    
    colors=zeros(nchoosek(K,3),3);
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                colors(trip_idx(i,j,k,K),:)=cell2mat(colors_c(i,j,k));
            end
        end
    end
    colors=sum(colors,2);
end


%criterion: product of ik*kj-ij is minimal
function [colors]=match_colors_c2(Rijs_rows,K)
    tic
    Rijs_rows_tpd=multi_transpose(Rijs_rows);
    trip_perms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
    inverse_perms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1];
    m=zeros(6,6);
    colors=zeros(nchoosek(K,3),3);
    votes=zeros(nchoosek(K,3),1);
    delta_arr=zeros(3,3,3,3,3);
    
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                prod_arr=abs(multiprod(permute(Rijs_rows(:,:,k3,:),[1,2,4,3]),Rijs_rows_tpd(:,:,k2,:)));
                delta_arr(:,:,1,:,:)=bsxfun(@minus,prod_arr,abs(Rijs_rows_tpd(:,:,k1,1)));
                delta_arr(:,:,2,:,:)=bsxfun(@minus,prod_arr,abs(Rijs_rows_tpd(:,:,k1,2)));
                delta_arr(:,:,3,:,:)=bsxfun(@minus,prod_arr,abs(Rijs_rows_tpd(:,:,k1,3)));
                norms=reshape(delta_arr,9,3,3,3);
                norms=squeeze(sum(norms.^2,1));
                for l=1:6
                    p1=trip_perms(l,:);
                    for r=1:6
                        p2=trip_perms(r,:);
                        m(l,r)=norms(1,p1(1),p2(1))+norms(2,p1(2),p2(2))+norms(3,p1(3),p2(3));
                    end
                end 
                idx=trip_idx(i,j,k,K);
                [votes(idx),tmp]=min(m(:));
                col=ceil(tmp/6);
                row=tmp-6*(col-1);
                colors(idx,1:2)=[row,col];
                
                %Need to calculate the relative permutation of Rik to Rij
                %by (sigma_ik)\circ(sigma_ij)^-1 
                inv_jk_perm=inverse_perms(col,:);
                rel_perm=inv_jk_perm(inv_jk_perm,:);
                colors(idx,3)=(2*rel_perm(1)-1)+(rel_perm(2)>rel_perm(3));
            end
        end
    end
    toc
end

function [sorted,classes]=sort_3(arr)
    idxs=zeros(3,1);
    sorted=zeros(3,1);
    
    if arr(1)<arr(2)
        idxs(1)=1;
        idxs(2)=2;
    else
        idxs(1)=2;%first member is second largest
        idxs(2)=1;
    end
    
    if arr(3)<arr(idxs(1))
        idxs(1:2)=idxs(1:2)+1;
        idxs(3)=1;
    elseif arr(3)<arr(idxs(2))
        idxs(idxs(2))=3;
        idxs(3)=2;
    else
        idxs(3)=3;
    end
    %sorted(idxs)=arr(idxs);
    [sorted,~]=sort(arr);
    classes=idxs;
        
end