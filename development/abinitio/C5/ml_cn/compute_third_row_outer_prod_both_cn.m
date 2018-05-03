function [vijs,viis,max_corrs_stats,mse_vij,detec_rate] = ...
    compute_third_row_outer_prod_both_cn(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,n_symm,max_shift,shift_step,is_handle_equators,refq,is_viz_cls)

if ~exist('is_viz_cls','var')
    is_viz_cls = false;
end

n_theta_ijs = size(R_theta_ijs,3);
n_theta_ijs_to_keep = floor(n_theta_ijs/n_symm)*n_symm;
if n_theta_ijs_to_keep < n_theta_ijs
    cijs_inds(:,:,n_theta_ijs_to_keep+1:end) = [];
    R_theta_ijs(:,:,n_theta_ijs_to_keep+1:end) = [];
    log_message('number of inplane rotation matrices must be divisible by n_symm=%d. Keeping %d/%d in-plane_rotation matrices',...
        n_symm,n_theta_ijs_to_keep,n_theta_ijs);
end


[n_r,n_theta,nImages] = size(npf);
nRis_tilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

log_message('precompiling shift phases');
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);

if ~exist('refq','var')
    refq = [];
end

log_message('precomputing likelihood with respect to self-common-lines');
all_self_corrs = estimate_viis(ciis,Ris_tilde,npf,n_symm,max_shift,shift_step,shift_phases,is_handle_equators,refq);
all_self_corrs(all_self_corrs == -inf) = 0;

log_message('computing likelihood with respect to all lines');
g_shift_phases = gpuArray(single(shift_phases));
[~,nshifts] = size(shift_phases);

opt_Rs_tilde    = zeros(nImages,nImages);
opt_thetaij     = zeros(nImages,nImages);
max_corrs_stats = zeros(nImages,nImages);


msg = [];

viis_tmp = zeros(3*nImages,3*nImages);
vijs = zeros(3*nImages,3*nImages);
g = [cosd(360/n_symm) -sind(360/n_symm) 0; 
	 sind(360/n_symm)  cosd(360/n_symm) 0; 
	 0 				 0  1]; % a rotation of 360/12 degrees about the z-axis
J = diag([1,1,-1]);
% g = gpuDevice(1);
% counter = 0;

% find which of the candidates rotations are equator images
equators_inds = find(abs(acosd(Ris_tilde(3,3,:)) - 90) < 7);
equators_inds = squeeze(equators_inds);

g_cijs_inds = gpuArray(single(cijs_inds));

for i=1:nImages
    % nomalize each ray to be norm 1
    npf_i = npf(:,:,i);
    norms   = sqrt(sum((abs(npf_i)).^2));
    npf_i = bsxfun(@rdivide,npf_i,norms);
    npf(:,:,i) = npf_i;
end

for i=1:nImages
    
    npf_i = npf(:,1:n_theta/2,i);
    g_npf_i = gpuArray(single(npf_i));
    
    g_npf_i_shifted = zeros([n_r,n_theta/2,nshifts],'gpuArray');
    for s=1:nshifts
        g_npf_i_shifted(:,:,s) = bsxfun(@times,g_npf_i,g_shift_phases(:,s));
    end
    
    g_npf_i_shifted = reshape(g_npf_i_shifted,n_r,n_theta/2*nshifts);
    
    
    g_all_self_corrs_i = gpuArray(single(all_self_corrs(:,i)));
    
    for j=i+1:nImages
        
        t1 = clock;
        
%         counter = counter + 1;
        
        
        g_all_self_corrs_j = gpuArray(single(all_self_corrs(:,j)));
        
        npf_j = npf(:,:,j);
        g_npf_j = gpuArray(single(npf_j));
        
        Corrs = g_npf_j'*g_npf_i_shifted; % Corrs is of dimension (n_theta x n_theta/2*nshifts )
        
        Corrs = reshape(Corrs,n_theta*n_theta/2,nshifts);
        Corrs = max(real(Corrs),[],2);
        Corrs = reshape(Corrs,n_theta,n_theta/2);
        Corrs = Corrs.';
        
        Corrs = Corrs(g_cijs_inds);

        Corrs = reshape(Corrs,[nRis_tilde,nRis_tilde,n_theta_ij/n_symm,n_symm]);
        
        %average over the common-lines
        Corrs = squeeze(mean(Corrs,4));
        
        [Corrs,KK] = max(Corrs,[],3);
        Corrs = Corrs.*(g_all_self_corrs_i*g_all_self_corrs_j.');
%         Corrs(Corrs < 0) = 0;
%         Corrs = log(Corrs.*(g_all_self_corrs_i*g_all_self_corrs_j.'));
%         Corrs = Corrs + bsxfun(@plus,g_all_self_corrs_i,g_all_self_corrs_j.');
        
        [YY,II] = max(Corrs(:));
        
        [II,JJ] = ind2sub([nRis_tilde,nRis_tilde],II);
        
        KK = KK(II,JJ);
        
        opt_Rs_tilde(i,j)    = gather(II);
        opt_Rs_tilde(j,i)    = gather(JJ);
        opt_thetaij(i,j)     = gather(KK);
        max_corrs_stats(i,j) = gather(YY);
        
        Ri_tilde = Ris_tilde(:,:,opt_Rs_tilde(i,j));
        Rj_tilde = Ris_tilde(:,:,opt_Rs_tilde(j,i));
        R_theta_ij = R_theta_ijs(:,:,opt_thetaij(i,j));
        
        vij = zeros(3,3);
        vii = zeros(3,3);
        vjj = zeros(3,3);
        for s=0:n_symm-1
            vij = vij + Ri_tilde.'*g^s*R_theta_ij*Rj_tilde;
            vii = vii + Ri_tilde.'*g^s*Ri_tilde;
            vjj = vjj + Rj_tilde.'*g^s*Rj_tilde;
        end        
        vijs((i-1)*3+1:i*3,(j-1)*3+1:j*3)     = vij/n_symm;
        viis_tmp((i-1)*3+1:i*3,(j-1)*3+1:j*3) = vii/n_symm;
        viis_tmp((j-1)*3+1:j*3,(i-1)*3+1:i*3) = vjj/n_symm;
        
        if true            
            %%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t2 = clock;
            t = etime(t2,t1);
            bs = char(repmat(8,1,numel(msg)));
            fprintf('%s',bs);
            msg = sprintf('i=%3d/%3d j=%3d/%3d t=%7.5f',i,nImages,j,nImages,t);
            fprintf('%s',msg);
            %%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
%     if mod(i,20) == 0
%         log_message('i=%3d/%3d',i,nImages);
%     end
end

viis = zeros(3,3,nImages);
for i=1:nImages
    
    vii_all = viis_tmp((i-1)*3+1:i*3,:);
    vii_all(:,(i-1)*3+1:i*3) = [];
    vii_all = reshape(vii_all,3,3,nImages-1);
    viis(:,:,i) = median(vii_all,3);
end

if false
    log_message('REMOVE THIS!!!!!!!! ONLY FOR DECUG PURPOSES');
    log_message('OVERRIDING VIJS');
    
    vijs = zeros(3,3,nchoosek(nImages,2));
    counter = 0;
    for i=1:nImages
        for j=i+1:nImages
            counter = counter + 1;
            
            rot_i_gt = q_to_rot(refq(:,i)).';
            rot_j_gt = q_to_rot(refq(:,j)).';
            
            vij = rot_i_gt(3,:).'* rot_j_gt(3,:);
            if rand < 0.5
                vij = J*vij*J;
            end
            vijs(:,:,counter) = vij;
        end
    end
end


if exist('refq','var') && ~isempty(refq)
    %% check mse of vij
    errs_mse_mat = zeros(nImages,nImages);
    errs_mse = zeros(1,nchoosek(nImages,2));
    counter = 0;
    for i=1:nImages
        for j=i+1:nImages
            counter = counter + 1;
            
            rot_i_gt = q_to_rot(refq(:,i)).';
            rot_j_gt = q_to_rot(refq(:,j)).';
            
            vij_gt = rot_i_gt(3,:).'* rot_j_gt(3,:);
            
            vij = vijs((i-1)*3+1:i*3,(j-1)*3+1:j*3);
            
            err = min(norm(vij-vij_gt,'fro'),norm(J*vij*J-vij_gt,'fro'));
            errs_mse_mat(i,j) = err;
            errs_mse(counter) = err;
        end
    end
    
    log_message('\n vij errs mse=%.2f',mean(errs_mse));
    
    
    
    %% check mse of vii
    errs_mse = zeros(1,nImages);
    for i=1:nImages
        
        rot_i_gt = q_to_rot(refq(:,i)).';
        
        vii_gt = rot_i_gt(3,:).'* rot_i_gt(3,:);
        
        vii = viis(:,:,i);
        
        err = min(norm(vii-vii_gt,'fro'),norm(J*vii*J-vii_gt,'fro'));
        errs_mse(i) = err;
        
    end
    
    log_message('\n vii errs mse=%.2f',mean(errs_mse));
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hand_idx = zeros(1,nchoosek(nImages,2));
    angle_tol_err = 10/180*pi; % how much angular deviation we allow for a common-line to have
    is_correct = zeros(nImages,nImages);
    counter = 0;
    for i=1:nImages
        
        for j=i+1:nImages
            
            counter = counter + 1;
            % step 1: compute the ground truth Ri*Rj and common-lines
            Ri_gt = q_to_rot(refq(:,i)).';
            Rj_gt = q_to_rot(refq(:,j)).';
            
            Rij_gt = Ri_gt.'*Rj_gt;
            
            c_1_gt = [-Rij_gt(8) ;  Rij_gt(7)]; %[-Riig(2,3)  Riig(1,3)];
            c_2_gt = [ Rij_gt(6) ; -Rij_gt(3)]; %[ Riig(3,2) -Riig(3,1)];
            
            cij_gt = clAngles2Ind(c_1_gt,n_theta);
            cji_gt = clAngles2Ind(c_2_gt,n_theta);
            
            
            % step 2: compute our estimate for Ri*g^sij*Rj
            Ri_tilde   = Ris_tilde(:,:,opt_Rs_tilde(i,j));
            Rj_tilde   = Ris_tilde(:,:,opt_Rs_tilde(j,i));
            R_theta_ij = R_theta_ijs(:,:,opt_thetaij(i,j));
            
            diff_s = zeros(1,n_symm);
            for s=1:n_symm
                Rij_g = Ri_tilde.'*g^(s-1)*R_theta_ij*Rj_tilde;
                    
                c_1 = [-Rij_g(8) ;  Rij_g(7)]; %[-Riig(2,3)  Riig(1,3)];
                c_2 = [ Rij_g(6) ; -Rij_g(3)]; %[ Riig(3,2) -Riig(3,1)];
                
                cij = clAngles2Ind(c_1,n_theta);
                cji = clAngles2Ind(c_2,n_theta);
                
                
                cij_diff = (cij_gt-cij)*2*pi./n_theta;
                cji_diff = (cji_gt-cji)*2*pi./n_theta;
                
                % take absolute cosine because of handedness
                % there might be +180 independendt diff for each image which at this stage
                % hasn't been taken care yet.
                diff = acos(cos(cij_diff))    + acos(cos(cji_diff));
                diff_J = acos(cos(cij_diff+pi)) + acos(cos(cji_diff+pi));
                if diff < diff_J
                    diff_s(s) = diff;
                    hand_idx(s) = 1;
                else
                    diff_s(s) = diff_J;
                    hand_idx(s) = 2;
                end
            end
            [min_diff_s,min_idx] = min(diff_s);
            if min_diff_s < 2*angle_tol_err
                is_correct(i,j) = 1;
                hand_idx(counter) = hand_idx(min_idx);
            end
        end
    end
    
    cl_dist = histc(hand_idx,1:2)/numel(hand_idx);
    detec_rate = sum(sum(is_correct))/nchoosek(nImages,2);
    log_message('\ncommon lines detection rate=%.2f%%',detec_rate*100);
    log_message('cl_J_dist=[%.2f %.2f]',cl_dist);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% find the angle between every two planes%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rel_rot_angles = zeros(nImages,nImages);
    for i=1:nImages
        for j=i+1:nImages
            Ri_gt = q_to_rot(refq(:,i)).';
            Rj_gt = q_to_rot(refq(:,j)).';
            
            angs = zeros(1,n_symm);
            for s=0:n_symm-1
                angs(s+1) = acosd(abs(Ri_gt(:,3).'*g^(s)*Rj_gt(:,3)));
            end
            rel_rot_angles(i,j) = min(angs); 
        end
    end
    % the bad indeces are those that common-lines were not found
    % correctly. put inf in lower part of matrix since these are all zeros
    bad_inds = find((is_correct + tril(inf*ones(nImages,nImages)))== 0);
    
%     figure; hist(rel_rot_angles(bad_inds)); title('histogram of angles between image planes whose common-lines were miss-detected');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualize those pair of images whose common-lines were not found and
    % the angle between the planes is not small%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:nImages
        for j=i+1:nImages
            if is_correct(i,j) == 0 && rel_rot_angles(i,j) > 10
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Ri_gt = q_to_rot(refq(:,i)).';
                Rj_gt = q_to_rot(refq(:,j)).';
                cijs_gt = zeros(1,n_symm);
                cjis_gt = zeros(1,n_symm);
                for s=1:n_symm
                    
                    RigsRj_gt = Ri_gt.'*g^(s-1)*Rj_gt;
                    
                    % extract a common-line index for each possible theta_ij
                    c_1 = [-RigsRj_gt(8) ;  RigsRj_gt(7)]; %[-U(2,3)  U(1,3)];
                    c_2 = [ RigsRj_gt(6) ; -RigsRj_gt(3)]; %[ U(3,2) -U(3,1)];
                    
                    c_1 = clAngles2Ind(c_1,n_theta);
                    c_2 = clAngles2Ind(c_2,n_theta);
                    
                    if c_2 > n_theta/2    
                        c_2 = c_2 - n_theta/2;
                        c_1 = c_1 + n_theta/2;
                        c_1 = mod(c_1-1,n_theta)+1;
                    end
                    
                    cijs_gt(s) = c_1;
                    cjis_gt(s) = c_2;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Ri_tilde   = Ris_tilde(:,:,opt_Rs_tilde(i,j));
                Rj_tilde   = Ris_tilde(:,:,opt_Rs_tilde(j,i));
                R_theta_ij = R_theta_ijs(:,:,opt_thetaij(i,j));
            
                cijs_found = zeros(1,n_symm);
                cjis_found = zeros(1,n_symm);
                for s=1:n_symm
                    
                    RigsRj = Ri_tilde.'*g^(s-1)*R_theta_ij*Rj_tilde;
                    
                    % extract a common-line index for each possible theta_ij
                    c_1 = [-RigsRj(8) ;  RigsRj(7)]; %[-U(2,3)  U(1,3)];
                    c_2 = [ RigsRj(6) ; -RigsRj(3)]; %[ U(3,2) -U(3,1)];
                    
                    c_1 = clAngles2Ind(c_1,n_theta);
                    c_2 = clAngles2Ind(c_2,n_theta);
                    
                    if c_2 > n_theta/2    
                        c_2 = c_2 - n_theta/2;
                        c_1 = c_1 + n_theta/2;
                        c_1 = mod(c_1-1,n_theta)+1;
                    end
                    
                    cijs_found(s) = c_1;
                    cjis_found(s) = c_2;
                end
                if is_viz_cls
                    viz_cls(cijs_found,cjis_found,cijs_gt,cjis_gt,n_theta);
                    pause;
                end
            end
        end
    end
    
    is_close_angle = rel_rot_angles + tril(inf*ones(size(rel_rot_angles)))< 10;
     
    both_eqtr_ims = false(nImages,nImages);
    im_eqtr_inds = [];
    for i=1:nImages
        Ri_gt = q_to_rot(refq(:,i)).';
        if(abs(acosd(Ri_gt(3,3)) - 90) < 7)
            im_eqtr_inds(end+1) = i;
        end
    end
    
    both_eqtr_ims(im_eqtr_inds,im_eqtr_inds) = true;
    is_both_eqtr_ims = both_eqtr_ims & tril(false*ones(size(both_eqtr_ims)));
    
    detec_rate_ex =  sum(sum(is_correct | is_close_angle | is_both_eqtr_ims))/nchoosek(nImages,2);
    log_message('\ncommon lines detection rate excluding close-planes and both planes are equators=%.2f%%',detec_rate_ex*100);
   
end

end