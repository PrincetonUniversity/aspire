function u_G = handedness_synchronization_TO(symmetry, est_rel_rots, cache_file_name)

% J sync method
%
% Parameters
%   symmetry          symmetry type: 'T' or 'O'.
%   est_rel_rots      array of estimated rotation matrices for each projection image.
% cache_file_name     The mat file name containing all candidate rotation
%                     matrices.  
%


% loading data
n_images = size(est_rel_rots,1);
gR = cryo_TO_group_elements(symmetry);
load(cache_file_name,'R');

n_gR = size(gR,3);

n_pairs = nchoosek(n_images,2);
Rijs = zeros(3,3,n_pairs,n_gR);
for i=1:n_images-1
    for j=i+1:n_images
        pair_ind = uppertri_ijtoind_vec(i,j,n_images);
        for k=1:n_gR
            Rijs(:,:,pair_ind,k) = R(:,:,est_rel_rots(i,j))*gR(:,:,k)*R(:,:,est_rel_rots(j,i)).';
        end
    end
end

u_G = jsync_TO(n_gR,Rijs,n_images);
end

% Synchronize handedness of relative rotations for T and O symmetrical molecules. 
function u_G = jsync_TO(n_gR,Rijs,K)
    log_message('Syncing relative rotations');
    J_list = outer_sync_brute_eff(n_gR,Rijs,K);
    log_message('Done matching triplets, now running power method...');  
    [u_G,~,evals] = signs_power_method (K,double(J_list),3,0); 
    log_message('top 3 eigenvalues = %.2f %.2f %.2f\n',evals(1),evals(2),evals(3));
    log_message('Done syncing handedness');
end

%
function [J_list]=outer_sync_brute_eff(n_gR,Rijs,K)

nworkers = 32;
J_list = zeros(nchoosek(K,3),1);
ntrip = nchoosek(K,3);
ntrip_per_worker = floor(ntrip/nworkers);
ntrip_last_worker = ntrip-ntrip_per_worker*(nworkers-1);
iter_lin_idx = zeros(nworkers+1,1);
for i=2:nworkers
    iter_lin_idx(i) = iter_lin_idx(i-1)+ntrip_per_worker;
end
iter_lin_idx(nworkers+1) = iter_lin_idx(nworkers)+ntrip_last_worker;

workers_res  = cell(1,nworkers);
parfor i=1:nworkers
    
    lin_idx_loc = iter_lin_idx;
    workers_res{i}=...
        outer_sync_brute_i(n_gR,Rijs,K,lin_idx_loc(i),lin_idx_loc(i+1));
end

for i=1:nworkers
    J_list(iter_lin_idx(i)+1:iter_lin_idx(i+1)) = workers_res{i};
end

end

%
function [J_list_i]=outer_sync_brute_i(n_gR,Rijs,K,from,to)

final_votes=zeros(4,1);
J_list_i=zeros(nchoosek(K,3),1);
norms=zeros(nchoosek(K,3),1);

ntrip=to-from;
trip_idx=lin2sub3_map(K);
trip_idx=trip_idx(from+1:to,:);

k1s=uppertri_ijtoind_vec(trip_idx(:,1),trip_idx(:,2),K);
k2s=uppertri_ijtoind_vec(trip_idx(:,1),trip_idx(:,3),K);
k3s=uppertri_ijtoind_vec(trip_idx(:,2),trip_idx(:,3),K);
ks=[k1s,k2s,k3s];
Rijs_t=permute(Rijs,[2,1,3,4]);

% For each triplets of indices i,j and k, consider the relative rotations
% tuples {Ri^TgmRj}, {Ri^TglRk} and {Rj^TgrRk}. Compute norms of the form
% ||Ri^TgmRj*Rj^TglRk-Ri^TglRk||, ||J*Ri^TgmRj*J*Rj^TglRk-Ri^TglRk||,
% ||Ri^TgmRj*J*Rj^TglRk*J-Ri^TglRk| and ||Ri^TgmRj*Rj^TglRk-J*Ri^TglRk*J||
% where gm,gl,gr are the varipus gorup members of Dn and J=diag([1,1-1]).
% The correct "J-configuration" is given for the samllest of these 4 norms.
for t=1:ntrip
    k1=ks(t,1);
    k2=ks(t,2);
    k3=ks(t,3);
    Rij=Rijs(:,:,k1,:);
    Rijk=cat(3,Rijs(:,:,k2,:),Rijs_t(:,:,k3,:));
    
    [final_votes(1),prod_arr]=compare_rot_brute_eff(n_gR,Rij,Rijk);
    final_votes(2)=compare_rot_brute_eff(n_gR,Rij,[],multi_Jify(prod_arr));
    k2_Jified=multi_Jify(Rijk(:,:,2,:));
    Rijk(:,:,2,:)=k2_Jified;
    [final_votes(3),prod_arr]=compare_rot_brute_eff(n_gR,Rij,Rijk);
    final_votes(4)=compare_rot_brute_eff(n_gR,Rij,[],multi_Jify(prod_arr));
    [norms(from+t),decision]=min(final_votes);
    J_list_i(from+t)=decision(1)-1;

end
J_list_i=J_list_i(from+1:to);

end

%
function [vote,prod_arr]=compare_rot_brute_eff(n_gR,Rij,Rijk,Jified_rot)
    %calculate distances between Rij estimates and original Rij
    if nargin<4
        prod_arr=multiprod(Rijk(:,:,1,:),squeeze(Rijk(:,:,2,:)));
    else
        prod_arr=Jified_rot;
    end
    
    arr=zeros(3,3,n_gR,n_gR,n_gR);
    for i=1:n_gR
        arr(:,:,:,:,i)=prod_arr-repmat(Rij(:,:,i),[1 1 n_gR n_gR]);
    end
    arr=reshape(arr,9,n_gR*n_gR*n_gR,1);
    arr=sum(arr.^2,1);
    [m,~]=sort(arr);
    vote=sum(m(1:n_gR*n_gR));
    
end