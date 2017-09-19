function [vijs,viis,max_corrs_stats,mse_vij,detec_rate] = ...
    compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,refq,is_viz_cls)
[n_r,n_theta,nImages] = size(npf);
nRis_tilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

%precompile the shift phases
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
viis = estimate_viis(ciis,Ris_tilde,npf,shift_phases,refq);

g_shift_phases = gpuArray(single(shift_phases));
[~,nshifts] = size(shift_phases);

opt_Rs_tilde    = zeros(nImages,nImages);
opt_thetaij     = zeros(nImages,nImages);
max_corrs_stats = zeros(nImages,nImages);

log_message('\ncomputing likelihood');
msg = [];

vijs = zeros(3,3,nchoosek(nImages,2));
g = [cosd(72) -sind(72) 0; 
	 sind(72)  cosd(72) 0; 
	 0 				 0  1]; % a rotation of 72 degrees about the z-axis
J = diag([1,1,-1]);
% g = gpuDevice(1);
counter = 0;

inds = sub2ind([n_theta,n_theta/2],cijs(:,:,:,:,1),cijs(:,:,:,:,2));
clear cijs;
[C,~,IC] = unique(inds(:));
clear inds;
g_C = gpuArray(single(C));
clear C;
g_IC = gpuArray(single(IC));
clear IC;
for i=1:nImages
    
    npf_i = npf(:,:,i);
    g_npf_i = gpuArray(single(npf_i));
    % ignoring dc-term
    g_npf_i(1,:) = 0;
    
    % nomalize each ray to be norm 1
    norms   = sqrt(sum((abs(g_npf_i)).^2));
    g_npf_i = bsxfun(@rdivide,g_npf_i,norms);
    
    for j=i+1:nImages
        
        t1 = clock;
        
        counter = counter + 1;
        
        npf_j = npf(:,1:n_theta/2,j);
        g_npf_j = gpuArray(single(npf_j));
        
        Corrs = zeros([numel(g_C),1],'gpuArray');
        for s=1:nshifts
            g_npf_j_shifted = bsxfun(@times,g_npf_j,g_shift_phases(:,s));
            
            % ignoring dc-term
            g_npf_j_shifted(1,:) = 0;
            
            % nomalize each ray to be norm 1
            norms           = sqrt(sum((abs(g_npf_j_shifted)).^2));
            g_npf_j_shifted = bsxfun(@rdivide,g_npf_j_shifted,norms);
            
            
            PiPj = g_npf_i'*g_npf_j_shifted;
            
            Corrs_s = PiPj(g_C);
            
            Corrs = max([Corrs real(Corrs_s(:))],[],2);
        end
        
        Corrs = real(Corrs(g_IC));
        Corrs = reshape(Corrs,[nRis_tilde,nRis_tilde,n_theta_ij/5,5]);
        %average over the five common-lines
        Corrs = squeeze(mean(Corrs,4));
        [YY,II] = max(Corrs(:));
        
        
        [II,JJ,KK] = ind2sub([nRis_tilde,nRis_tilde,n_theta_ij/5],II);
        
        opt_Rs_tilde(i,j)    = gather(II);
        opt_Rs_tilde(j,i)    = gather(JJ);
        opt_thetaij(i,j)     = gather(KK);
        max_corrs_stats(i,j) = gather(YY);
        
        Ri_tilde = Ris_tilde(:,:,opt_Rs_tilde(i,j));
        Rj_tilde = Ris_tilde(:,:,opt_Rs_tilde(j,i));
        R_theta_ij = R_theta_ijs(:,:,opt_thetaij(i,j));
        
        vij = zeros(3,3);
        for s=0:4
            vij = vij + Ri_tilde.'*g^s*R_theta_ij*Rj_tilde;
        end
        vijs(:,:,counter) = vij/5;
        
        %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t2 = clock;
        t = etime(t2,t1);
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('i=%3d/%3d j=%3d/%3d t=%7.5f',i,nImages,j,nImages,t);
        fprintf('%s',msg);
        %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

if exist('refq','var') && ~isempty(refq)
    
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
            
            diff_s = zeros(1,5);
            for s=1:5
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
            
            rel_rot_angles(i,j) = min([acosd(abs(Ri_gt(:,3).'*Rj_gt(:,3))),...
                acosd(abs(Ri_gt(:,3).'*g*Rj_gt(:,3))),...
                acosd(abs(Ri_gt(:,3).'*g^2*Rj_gt(:,3))),...
                acosd(abs(Ri_gt(:,3).'*g^3*Rj_gt(:,3)))...
                acosd(abs(Ri_gt(:,3).'*g^4*Rj_gt(:,3)))]);
            
        end
    end
    % the bad indeces are those that common-lines were not found
    % correctly. put inf in lower part of matrix since these are all zeros
    bad_inds = find((is_correct + tril(inf*ones(nImages,nImages)))== 0);
    
    figure; hist(rel_rot_angles(bad_inds)); title('histogram of angles between image planes whose common-lines were miss-detected');
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
                cijs_gt = zeros(1,5);
                cjis_gt = zeros(1,5);
                for s=1:5
                    
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
            
                cijs_found = zeros(1,5);
                cjis_found = zeros(1,5);
                for s=1:5
                    
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
    detec_rate_ex_close_planes =  sum(sum(or(is_correct,is_close_angle)))/nchoosek(nImages,2);
    log_message('\ncommon lines detection rate excluding close planes =%.2f%%',detec_rate_ex_close_planes*100);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the mse of outer-product of third rows
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    diffs = zeros(1,nchoosek(nImages,2));
    for i=1:nImages
        
        for j=i+1:nImages
            
            Ri_gt = q_to_rot(refq(:,i)).';
            Rj_gt = q_to_rot(refq(:,j)).';
            Rij_gt = Ri_gt.'*Rj_gt;
            % step 2: compute our estimate for Ri*g^sij*Rj
            Ri_tilde   = Ris_tilde(:,:,opt_Rs_tilde(i,j));
            Rj_tilde   = Ris_tilde(:,:,opt_Rs_tilde(j,i));
            R_theta_ij = R_theta_ijs(:,:,opt_thetaij(i,j));
            
            diff = zeros(1,5);
            for s=0:4
                Rij_g = Ri_tilde.'*g^s*R_theta_ij*Rj_tilde;
                diff(s+1) = min(norm(Rij_g-Rij_gt,'fro'),norm(J*Rij_g*J-Rij_gt,'fro'));
            end
            diffs(counter) = min(diff);
        end
    end
    
    mse_vij = sum(diffs.^2)/numel(diffs);
    log_message('MSE of vij: %e',mse_vij);
end

end