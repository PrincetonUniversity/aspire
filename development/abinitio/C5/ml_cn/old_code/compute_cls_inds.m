function [cijs_inds,Ris_tilde,R_theta_ijs,n_theta] = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name)

nRisTilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

cijs = zeros(nRisTilde,nRisTilde,n_theta_ij,2,'uint16');
msg = [];
for i=1:nRisTilde
    Ri_tilde = Ris_tilde(:,:,i);
    % calculate R_i_tilde*R_theta_ij for all possible R_theta_ijs
    tmp = multiprod(Ri_tilde.',R_theta_ijs);
    for j=1:nRisTilde
        t1 = clock;
        Rj_tilde = Ris_tilde(:,:,j);
        Rij = multiprod(tmp,Rj_tilde);
        
        Rij = reshape(Rij,9,n_theta_ij);
        % extract a common-line index for each possible theta_ij
        c_1 = [-Rij(8,:) ;  Rij(7,:)]; %[-U(2,3)  U(1,3)];
        c_2 = [ Rij(6,:) ; -Rij(3,:)]; %[ U(3,2) -U(3,1)];
        
        c_1s = clAngles2Ind(c_1,n_theta);
        c_2s = clAngles2Ind(c_2,n_theta);
        
        inds_tmp = find(c_1s > n_theta/2);
        c_1s(inds_tmp) = c_1s(inds_tmp) - n_theta/2;
        c_2s(inds_tmp) = c_2s(inds_tmp) + n_theta/2;
        c_2s(inds_tmp) = mod(c_2s(inds_tmp)-1,n_theta)+1;
        
        cijs(i,j,:,1) = c_1s;
        cijs(i,j,:,2) = c_2s;
        
        %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t2 = clock;
        t = etime(t2,t1);
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('i=%3d/%3d j=%3d/%3d t=%7.5f',i,nRisTilde,j,nRisTilde,t);
        fprintf('%s',msg);
        %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

cijs_inds = uint16(sub2ind([n_theta/2,n_theta],cijs(:,:,:,1),cijs(:,:,:,2)));
clear cijs;

if is_save_inds_to_cache
    log_message('saving inds to cache file=%s',file_cache_name);
    save(file_cache_name,'cijs_inds','Ris_tilde','R_theta_ijs','n_theta','-v7.3');
end
end