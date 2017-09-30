function [ciis,cijs] = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name)

nRisTilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

ciis = zeros(2,nRisTilde,'uint16');
g = [0 -1 0; 1 0 0; 0 0 1]; % a rotation of 90 degrees about the z-axis
for i=1:nRisTilde
    
    Ri_tilde = Ris_tilde(:,:,i);
    Riig = Ri_tilde.' * g * Ri_tilde;
    
    % extract a common-line index for each possible theta_ij
    c_1 = [-Riig(8) ;  Riig(7)]; %[-Riig(2,3)  Riig(1,3)];
    c_2 = [ Riig(6) ; -Riig(3)]; %[ Riig(3,2) -Riig(3,1)];
    
    cii1 = clAngles2Ind(c_1,n_theta);
    cii2 = clAngles2Ind(c_2,n_theta);
    
    if cii2 > n_theta/2
        cii2 = cii2 - n_theta/2;
        cii1 = cii1 + n_theta/2;
        cii1 = mod(cii1-1,n_theta)+1;
    end
    
    ciis(1,i)  = cii1;
    ciis(2,i)  = cii2;
end

cijs = zeros(nRisTilde,nRisTilde,n_theta_ij/4,4,2,'uint16');
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
        
        inds_tmp = find(c_2s > n_theta/2);
        c_2s(inds_tmp) = c_2s(inds_tmp) - n_theta/2;
        c_1s(inds_tmp) = c_1s(inds_tmp) + n_theta/2;
        c_1s(inds_tmp) = mod(c_1s(inds_tmp)-1,n_theta)+1;
        
        c_1s = reshape(c_1s,n_theta_ij/4,4);
        c_2s = reshape(c_2s,n_theta_ij/4,4);
        
        cijs(i,j,:,:,1) = c_1s;
        cijs(i,j,:,:,2) = c_2s;
        
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
if is_save_inds_to_cache
    save(file_cache_name,'ciis','cijs','Ris_tilde','R_theta_ijs','-v7.3');
end

end
