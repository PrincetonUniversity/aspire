
function [Rijs_synced,J_list,evals]=outer_sync_dev(Rijs,K,scheme,J_list_in)
    disp('outer syncing');
    tic
    %Rijs=reshape(in,3,3,nchoosek(K,2),4);
    third_rows=reshape(Rijs(3,:,:,:),3,nchoosek(K,2),4);
    
    switch scheme
        case 1
        %unweighted voting
        [J_list]=outer_sync_unweighted(Rijs,third_rows,K);
        case 2
        %weighted voting
        [J_list]=outer_sync_weighted1(Rijs,third_rows,K);
        case 3
        %weighted with adjusted voting
        [J_list]=outer_sync_weighted2(Rijs,third_rows,K);
        case 4
        %plain brute force
        [J_list]=outer_sync_brute1(Rijs,K);   
        case 5
        %smart brute force?
        [J_list]=outer_sync_brute2(Rijs,K);
        case 6
        [J_list,norms]=parallel_outer_sync_brute1(Rijs,K);
        case 7
        %Efficient version of parallel_outer_sync_brute1 
        [J_list]=outer_sync_brute_eff(Rijs,K,12);
        case 8
        J_list=J_list_in;
    end
    Rijs_synced=Rijs;
    toc
    disp('finished matching triplets');
%     saveDir='/a/home/cc/math/eitanros/cryo/Pipeline/Results';
%     save([saveDir,'/J_list_',cell2mat(datetime_suffix)],'J_list','-v7.3');
    tic
%     %FOR DEBUG
%     d=nchoosek(K,2);
%     sync_mat=zeros(d,d);
%     for i=1:K-2
%         for j=i+1:K-1
%             for k=j+1:K
%                 conf=J_list(trip_idx(i,j,k,K));
%                 k1=uppertri_ijtoind(i,j,K);
%                 k2=uppertri_ijtoind(i,k,K);
%                 k3=uppertri_ijtoind(j,k,K);
%                 sync_mat([k1,k2,k3],[k1,k2,k3])=ones(3,3);
%                 switch conf
%                     case 1
%                         sync_mat(k1,[k2,k3])=-1;
%                         sync_mat([k2,k3],k1)=-1;
%                     case 2
%                         sync_mat(k2,[k1,k3])=-1;
%                         sync_mat([k1,k3],k2)=-1;
%                     case 3
%                         sync_mat(k3,[k1,k2])=-1;
%                         sync_mat([k1,k2],k3)=-1;
%                 end                    
%             end
%         end
%     end
%     sync_mat=sync_mat-diag(ones(d,1));
%     [e_vec,e_val]=eigs(sync_mat,1);
%      Rijs_synced(:,:,e_vec<0,:)=multi_Jify(Rijs_synced(:,:,e_vec<0,:));
%      [J_list_new]=outer_sync_brute1(Rijs_synced,50);  
%      
%      %END DEBUG
    %save('J_list.mat','J_list');
    [J_signs,~,evals] =signs_power_method (K,double(J_list),3,0);
    %save('J_signs.mat','J_signs');
    %shift J_signs
%     neg_avg=average(J_signs(J_signs<0));
%     pos_avg=average(J_signs(J_signs>0));
%    J_signs=J_signs-average(J_signs);    
    
    Rijs_synced(:,:,J_signs<0,:)=multi_Jify(Rijs_synced(:,:,J_signs<0,:));
    
%     m=zeros(nchoosek(K,2),1);
%     for i=1:K-2
%         for j=i+1:K-1
%             for k=j+1:K
%                 ij=uppertri_ijtoind(i,j,K);
%                 ik=uppertri_ijtoind(i,k,K);
%                 jk=uppertri_ijtoind(j,k,K);
%                 [m(ij)]=pair_variance(Rijs_synced(:,:,[ik,jk],:));
%             end
%         end
%     end
    
    toc
    %[J_list_new]=outer_sync_brute1(Rijs_synced,K);  
end

function [J_list]=outer_sync_brute_eff(Rijs,K,nworkers)

J_list=zeros(nchoosek(K,3),1);
ntrip=nchoosek(K,3);
ntrip_per_worker=floor(ntrip/nworkers);
ntrip_last_worker=ntrip-ntrip_per_worker*(nworkers-1);
iter_lin_idx=zeros(nworkers+1,1);
for i=2:nworkers
    iter_lin_idx(i)=iter_lin_idx(i-1)+ntrip_per_worker;
end
iter_lin_idx(nworkers+1)=iter_lin_idx(nworkers)+ntrip_last_worker;

workers_res=cell(1,nworkers);
%h=waitbar(0,'outer syncing');

parfor i=1:nworkers
    
    lin_idx_loc=iter_lin_idx;
    workers_res{i}=...
        outer_sync_brute_i(Rijs,K,lin_idx_loc(i),lin_idx_loc(i+1));    
    %waitbar(i/nworkers);
end
%close(h);

for i=1:nworkers
    J_list(iter_lin_idx(i)+1:iter_lin_idx(i+1))=workers_res{i};
end

end

function [J_list_i]=outer_sync_brute_i(Rijs,K,from,to)

final_votes=zeros(4,1);
J_list_i=zeros(nchoosek(K,3),1);
norms=zeros(nchoosek(K,3),1);

ntrip=to-from;
trip_idx=lin2sub3_map(K);
trip_idx=trip_idx(from+1:to,:);

k1s=uppertri_ijtoind(trip_idx(:,1),trip_idx(:,2),K);
k2s=uppertri_ijtoind(trip_idx(:,1),trip_idx(:,3),K);
k3s=uppertri_ijtoind(trip_idx(:,2),trip_idx(:,3),K);
ks=[k1s,k2s,k3s];
Rijs_t=permute(Rijs,[2,1,3,4]);

for t=1:ntrip
%     [i,j,k]=idx_map(t,:);
%     k1=uppertri_ijtoind(i,j,K);
%     k2=uppertri_ijtoind(i,k,K);
%     k3=uppertri_ijtoind(j,k,K);
    k1=ks(t,1);
    k2=ks(t,2);
    k3=ks(t,3);
    Rij=Rijs(:,:,k1,:);
    %Rijk=Rijs(:,:,[k2 k3],:);
    Rijk=cat(3,Rijs(:,:,k2,:),Rijs_t(:,:,k3,:));
    
    [final_votes(1),prod_arr]=compare_rot_brute_eff(Rij,Rijk);
    final_votes(2)=compare_rot_brute_eff(Rij,[],multi_Jify(prod_arr));
    k2_Jified=multi_Jify(Rijk(:,:,2,:));
    Rijk(:,:,2,:)=k2_Jified;
    [final_votes(3),prod_arr]=compare_rot_brute_eff(Rij,Rijk);
    final_votes(4)=compare_rot_brute_eff(Rij,[],multi_Jify(prod_arr));
    [norms(from+t),decision]=min(final_votes);
    J_list_i(from+t)=decision(1)-1;

end
J_list_i=J_list_i(from+1:to);

end

%Plain brute force
function [vote,prod_arr]=compare_rot_brute_eff(Rij,Rijk,Jified_rot)
    %calculate distances between Rij estimates and original Rij
    if nargin<3
        prod_arr=multiprod(Rijk(:,:,1,:),squeeze(Rijk(:,:,2,:)));
    else
        prod_arr=Jified_rot;
    end
    
    arr=zeros(3,3,8,8);
    arr(:,:,1:4,1:4)=prod_arr-repmat(Rij(:,:,1),[1 1 4 4]);
    arr(:,:,1:4,5:8)=prod_arr-repmat(Rij(:,:,2),[1 1 4 4]);
    arr(:,:,5:8,1:4)=prod_arr-repmat(Rij(:,:,3),[1 1 4 4]);
    arr(:,:,5:8,5:8)=prod_arr-repmat(Rij(:,:,4),[1 1 4 4]);
    arr=reshape(arr,9,64,1);
    arr=sum(arr.^2,1);
    [m,~]=sort(arr);
    vote=sum(m(1:16));
    
end

function [J_list]=outer_sync_brute1(Rijs,K)
    tic
    final_votes=zeros(4,1);
    J_list=zeros(nchoosek(K,3),1);
    dummy=0;
    norms=zeros(nchoosek(K,3),1);
    h=waitbar(0,'J-syncing...');
    ijk_idx=0;
    ntrip=nchoosek(K,3);
    
    for i=1:K-2
        waitbar(ijk_idx/ntrip);
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                Rij=Rijs(:,:,k1,:);
                Rijk=Rijs(:,:,[k2 k3],:);       
                
                [final_votes(1),prod_arr]=compare_rot_brute1(Rij,Rijk);                   
                final_votes(2)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,2,:));
                Rijk(:,:,2,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_brute1(Rij,Rijk);           
                final_votes(4)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr)); 
                [norms(trip_idx(i,j,k,K)),decision]=min(final_votes);              
                J_list(trip_idx(i,j,k,K))=decision(1)-1;
                ijk_idx=ijk_idx+1;
            end
        end
    end
    toc
    close(h);
end

function [out,norms]=parallel_outer_sync_brute1(Rijs,K)
    tic
    
    J_list=cell(K,K,K); 
    norms_cell=cell(K,K,K);
    norms=zeros(nchoosek(K,3),1);
    
    parfor i=1:K-2
        in=Rijs;
        final_votes=zeros(4,1);
        dummy=0;
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
                Rij=in(:,:,k1,:);
                Rijk=in(:,:,[k2 k3],:);  
                %[j_val,~]=brute_iter_ijk(Rij,Rijk);
                %J_list(i,j,k)={int8(j_val)};
                
                [final_votes(1),prod_arr]=compare_rot_brute1(Rij,Rijk);
                final_votes(2)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,2,:));
                Rijk(:,:,2,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_brute1(Rij,Rijk);
                final_votes(4)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr));
                [tmp,decision]=min(final_votes);
                J_list(i,j,k)={decision(1)-1};
                norms_cell(i,j,k)={tmp};
            end
        end
    end
    out=zeros(nchoosek(K,3),1);
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                out(trip_idx(i,j,k,K))=cell2mat(J_list(i,j,k));
                norms(trip_idx(i,j,k,K))=cell2mat(norms_cell(i,j,k));
            end
        end
    end
    toc
end


%Plain brute force
function [vote,prod_arr]=compare_rot_brute1(Rij,Rijk,Jified_rot)
    %calculate distances between Rij estimates and original Rij
    if nargin<3
        trans_rot3=reshape(Rijk(:,:,2,:),9,4);
        trans_rot3=trans_rot3([1 4 7 2 5 8 3 6 9],:);
        trans_rot3=reshape(trans_rot3,3,3,4);
        %rot3=reshape(rot3,3,3,1,4);%TO DO: Check if needed
        prod_arr=multiprod(Rijk(:,:,1,:),trans_rot3);
    else
        prod_arr=Jified_rot;
    end
    
    arr=zeros(3,3,8,8);
    arr(:,:,1:4,1:4)=prod_arr-repmat(Rij(:,:,1),[1 1 4 4]);
    arr(:,:,1:4,5:8)=prod_arr-repmat(Rij(:,:,2),[1 1 4 4]);
    arr(:,:,5:8,1:4)=prod_arr-repmat(Rij(:,:,3),[1 1 4 4]);
    arr(:,:,5:8,5:8)=prod_arr-repmat(Rij(:,:,4),[1 1 4 4]);
    arr=reshape(arr,9,64,1);
    arr=sum(arr.^2,1);
    [m,~]=sort(arr);
    vote=sum(m(1:16));
    
    %TO DO: Try best feasible option. 
    
    %Frobenius multi-norm
%     prod_arr=reshape(prod_arr,9,4,4);
%     prod_arr=prod_arr.^2;
%     m=sum(prod_arr,1);
%     vote=sum(m);
    
end

function [j_val,norms]=brute_iter_ijk(Rij,Rijk)
    dummy=0;
    final_votes=zeros(4,1);
    [final_votes(1),prod_arr]=compare_rot_brute1(Rij,Rijk);
    final_votes(2)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr));
    k2_Jified=multi_Jify(Rijk(:,:,2,:));
    Rijk(:,:,2,:)=k2_Jified;
    [final_votes(3),prod_arr]=compare_rot_brute1(Rij,Rijk);
    final_votes(4)=compare_rot_brute1(Rij,dummy,multi_Jify(prod_arr));
    [norms,decision]=min(final_votes);
    j_val=decision(1)-1;
end


function [J_list]=outer_sync_brute2(Rijs,K)
    tic
    final_votes=zeros(4,1);
    J_list=zeros(nchoosek(K,3),1);
    p=perms([1 2 3 4]);
    perm_Rijs=zeros([size(Rijs) 24]);
    mega_perm=1:96;
    mega_perm=reshape(mega_perm,24,4);
    mega_perm=reshape(mega_perm',96,1);
    norms=zeros(nchoosek(K,3),1);%FOR DEBUG

    for i=1:24
        perm_Rijs(:,:,:,:,i)=Rijs(:,:,:,p(i,:));
    end
    
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K

                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                Rijk=Rijs(:,:,[k2 k3],:);  
                Rij=squeeze(perm_Rijs(:,:,k1,:,:));
                
                [final_votes(1),prod_arr]=compare_rot_brute2(Rij,Rijk,mega_perm);                   
                final_votes(2)=compare_rot_brute2(Rij,Rijk,mega_perm,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,1,:));
                Rijk(:,:,1,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_brute2(Rij,Rijk,mega_perm);           
                final_votes(4)=compare_rot_brute2(Rij,Rijk,mega_perm,multi_Jify(prod_arr)); 
                [norms(trip_idx(i,j,k,K)),desicion]=min(final_votes);                                
                J_list(trip_idx(i,j,k,K))=desicion(1)-1;
                
            end
        end
    end
    toc
end

function [vote,prod_arr]=compare_rot_brute2(Rij,Rijk,p,Jified_rot)
    %calculate distances between Rij estimates and original Rij
    if nargin<4
        trans_rot3=reshape(Rijk(:,:,2,:),9,4);
        trans_rot3=trans_rot3([1 4 7 2 5 8 3 6 9],:);
        trans_rot3=reshape(trans_rot3,3,3,4);
        %rot3=reshape(rot3,3,3,1,4);%TO DO: Check if needed
        prod_arr=multiprod(Rijk(:,:,1,:),trans_rot3);
        
        prod_arr=permute(prod_arr,[4 3 1 2]);
        prod_arr=repmat(prod_arr,24,1);
        prod_arr=permute(prod_arr,[3 4 2 1]);
    else
        prod_arr=Jified_rot;
    end
    
    %multi repmat inline
    Rijs_rep=permute(Rij,[4 3 1 2]);
    Rijs_rep=repmat(Rijs_rep,4,1);
    Rijs_rep=permute(Rijs_rep,[3 4 2 1]);
    Rijs_rep=Rijs_rep(:,:,:,p);
    
    diff=prod_arr-Rijs_rep;
    diff=reshape(diff,9,4,4,24);
    diff=sum(diff.^2,1);%Frobenius multi-norm
     
    aggregate_votes=squeeze(sum(diff,2));
    %aggregate_votes=reshape(aggregate_votes,4,24);
    [agg_weights,agg_votes]=min(aggregate_votes);
    
    %for each of 4 Rijk options find which of 24 permutations of Rij is 
    %closest by norm
    perm_idxs=1:24;
    voters=unique(agg_votes);
    best_perms=zeros(length(voters),1);
    count=0;
    for i=voters
        count=count+1;
        i_voters=agg_votes==i;%Can be Empty
        [~,idx]=min(agg_weights(i_voters));
        i_voters=perm_idxs(i_voters);
        best_perms(count)=i_voters(idx);
    end
%     if ~isinteger(best_perms)
%         vote=0;
%     end
    vote=sum(agg_weights(best_perms));
    
%     %weight the votes with adjusted weights
%     adj_weights=weights_list;
%     for i=1:4
%         adj_weights(i,:)=agg_weights(i)*adj_weights(i,:);
%     end
%     m=weights_list.*m;
%     vote=sum(sum(m)); 
end

function [J_list]=outer_sync_unweighted(Rijs,third_rows,K)
            
    final_votes=zeros(4,1);
    J_list=zeros(nchoosek(K,3),1);
    p=perms(1:4);
    
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                Rij=Rijs(:,:,k1,:);
                Rijk=Rijs(:,:,[k2 k3],:);       
                [votes_list1,sum_weights1]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                [votes_list2,sum_weights2]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                if sum_weights1>sum_weights2
                    votes_list=votes_list1;
                else
                    votes_list=votes_list2;
                end
                
                [final_votes(1),prod_arr]=compare_rot_unweighted(votes_list,Rij,Rijk);                   
                final_votes(2)=compare_rot_unweighted(votes_list,Rij,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,1,:));
                Rijk(:,:,1,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_unweighted(votes_list,Rij,Rijk);           
                final_votes(4)=compare_rot_unweighted(votes_list,Rij,multi_Jify(prod_arr)); 
                [~,decision]=min(final_votes);                                
                J_list(trip_idx(k1,k2,k3,K))=desicion(1)-1;
                
            end
        end
    end
end

function [J_list]=outer_sync_weighted1(Rijs,third_rows,K)
            
    final_votes=zeros(4,1);
    J_list=zeros(nchoosek(K,3),1);
    p=perms(1:4);
    
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                Rij=Rijs(:,:,k1,:);
                Rijk=Rijs(:,:,[k2 k3],:);       
                [votes_list1,weights_list1]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                [votes_list2,weights_list2,sum_weights2]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                if sum_weights1>sum_weights2
                    votes_list=votes_list1;
                    weights_list=weights_list1;
                else
                    votes_list=votes_list2;
                    weights_list=weights_list2;
                end
                
                [final_votes(1),prod_arr]=compare_rot_weighted1(votes_list,weights_list,Rij,Rijk);                   
                final_votes(2)=compare_rot_weighted1(votes_list,weights_list,Rij,0,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,1,:));
                Rijk(:,:,1,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_weighted1(votes_list,weights_list,Rij,Rijk);           
                final_votes(4)=compare_rot_weighted1(votes_list,weights_list,Rij,0,multi_Jify(prod_arr)); 
                [~,desicion]=min(final_votes);                                
                J_list(trip_idx(k1,k2,k3,K))=desicion(1)-1;
                
            end
        end
    end
end

function [J_list]=outer_sync_weighted2(Rijs,third_rows,K)
            
    final_votes=zeros(4,1);
    J_list=zeros(nchoosek(K,3),1);
    p=perms(1:4);
    
    for i=1:K-2
        for j=i+1:K-1
            for k=j+1:K
                k1=uppertri_ijtoind(i,j,K);
                k2=uppertri_ijtoind(i,k,K);
                k3=uppertri_ijtoind(j,k,K);
                Rij=Rijs(:,:,k1,:);
                Rijk=Rijs(:,:,[k2 k3],:);       
                [votes_list1,weights_list1]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                [votes_list2,weights_list2,sum_weights2]=best_perm_dep(third_rows(3,k1,:),third_rows(:,[k2,k3],:),p);
                if sum_weights1>sum_weights2
                    votes_list=votes_list1;
                    weights_list=weights_list1;
                else
                    votes_list=votes_list2;
                    weights_list=weights_list2;
                end
                
                [final_votes(1),prod_arr]=compare_rot_weighted2(votes_list,weights_list,Rij,Rijk);                   
                final_votes(2)=compare_rot_weighted2(votes_list,weights_list,agg_weights,Rij,0,multi_Jify(prod_arr));
                k2_Jified=multi_Jify(Rijk(:,:,1,:));
                Rijk(:,:,1,:)=k2_Jified;
                [final_votes(3),prod_arr]=compare_rot_weighted2(votes_list,weights_list,Rij,Rijk);           
                final_votes(4)=compare_rot_weighted2(votes_list,weights_list,agg_weights,Rij,0,multi_Jify(prod_arr)); 
                [~,desicion]=min(final_votes);                                
                J_list(trip_idx(k1,k2,k3,K))=desicion(1)-1;
                
            end
        end
    end
end

%Choose 16 best options out of 64. Minimizes 3,3 entry of |Rij-Rik*Rkj|.
function [votes_list,sum_weights,weights_list,agg_weights]=best_perm_dep(ij_angles,third_rows,p)

    %calculate 16 relative rotation angles ( entry 3,3 of Rij)
    prods=third_rows(:,1,:)*third_rows(:,2,:)';
       
    %24 different permutations of the 4 Rij
    diff=zeros(4,24*4);
    for i=1:24
        diff(1:4,4*i-3:4*i)=abs(prods-ij_angles(p(i,:)));
    end
    aggregate_votes=sum(diff,2);
    aggregate_votes=reshape(aggregate_votes,4,24);
    [agg_weights,agg_votes]=min(aggregate_votes,1);
    
    %create list of voters for all Rij
    %votes_list is 4x4 array where cell l,m contains the index of the Rij
    %which is closest to (l'th Rik)*(m'th Rkj)
    agg_votes=agg_votes(1,:);
    votes_list=p(agg_votes,:);%TO DO: maybe transpose
    weights_list=diff(agg_votes,:);   
    sum_weights=sum(agg_votes);
end

function [vote,prod_arr]=compare_rot_unweighted(votes_list,Rij,Rijk,Jified_rot)

    %calculate distances between Rij estimates and original Rij
    if nargin<4
        trans_rot3=reshape(Rijk(:,:,2,:),9,4);
        trans_rot3=trans_rot3([1 4 7 2 5 8 3 6 9],:);
        trans_rot3=reshape(trans_rot3,3,3,4);
        %rot3=reshape(rot3,3,3,1,4);%TO DO: Check if needed
        prod_arr=multiprod(reshape(Rijk(:,:,1,:),3,3,1,4),trans_rot3);
    else
        prod_arr=Jified_rot;
    end
    
    for i=1:4
        prod_arr(:,:,i,:)=prod_arr(:,:,i,:)-Rij(:,:,1,votes_list(i,:));
    end
    %Frobenius multi-norm
    prod_arr=reshape(prod_arr,9,4,4);
    prod_arr=prod_arr.^2;
    m=sum(prod_arr,1);
    vote=sum(m);
end

function [vote]=compare_rot_weighted1(votes_list,Rij,Rijk,weights_list,Jified_rot)

    %calculate distances between Rij estimates and original Rij
    if nargin<5
        trans_rot3=reshape(Rijk(:,:,2,:),9,4);
        trans_rot3=trans_rot3([1 4 7 2 5 8 3 6 9],:);
        trans_rot3=reshape(trans_rot3,3,3,4);
        %rot3=reshape(rot3,3,3,1,4);%TO DO: Check if needed
        prod_arr=multiprod(reshape(Rijk(:,:,1,:),3,3,1,4),trans_rot3);
    else
        prod_arr=Jified_rot;
    end
    
    for i=1:4
        prod_arr(:,:,i,:)=prod_arr(:,:,i,:)-Rij(:,:,votes_list(i,:));
    end
    %Frobenius multi-norm
    prod_arr=reshape(prod_arr,9,4,4);
    prod_arr=prod_arr.^2;
    m=sum(prod_arr,1);
    
    %weight the votes
    m=weights_list.*m;
    vote=sum(sum(m));
    
end

function [vote]=compare_rot_weighted2(votes_list,Rij,Rijk,weights_list,agg_weights,Jified_rot)

    %calculate distances between Rij estimates and original Rij
    if nargin<6
        trans_rot3=reshape(Rijk(:,:,2,:),9,4);
        trans_rot3=trans_rot3([1 4 7 2 5 8 3 6 9],:);
        trans_rot3=reshape(trans_rot3,3,3,4);
        %rot3=reshape(rot3,3,3,1,4);%TO DO: Check if needed
        prod_arr=multiprod(reshape(Rijk(:,:,1,:),3,3,1,4),trans_rot3);
    else
        prod_arr=Jified_rot;
    end
    
    for i=1:4
        prod_arr(:,:,i,:)=prod_arr(:,:,i,:)-Rij(:,:,votes_list(i,:));
    end
    %Frobenius multi-norm
    prod_arr=reshape(prod_arr,9,4,4);
    prod_arr=prod_arr.^2;
    m=sum(prod_arr,1);
    
    %weight the votes with adjusted weights
    adj_weights=weights_list;
    for i=1:4
        adj_weights(i,:)=agg_weights(i)*adj_weights(i,:);
    end
    m=weights_list.*m;
    vote=sum(sum(m));    
end

%Jify all matrices in multi-dim array assuming first 2 dimensions are the
%matrices
% function [Jified_rot]=multi_Jify(in)
%     
%     dims=size(in);
%     l=len(dims);
%     Jified_rot=reshape(in,dims(1)*dims(2),dims(3:l));
%     Jified_rot([3 6 7 8],:,:)=-Jified_rot([3 6 7 8],dims(3:l));
%     Jified_rot=reshape(Jified_rot,3,3,dims(3:l));
% 
% end

function [arr_res]=multi_repmat(arr)

    arr_p=permute(arr,[4 3 1 2]);
    arr_rep=repmat(arr_p,1,4);
    arr_res=permute(arr_rep,[3 4 1 2]);

end
