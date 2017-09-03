function [vijs,viis,max_corrs_stats] = ...
    compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,R_theta_ijs,refq)

[~,n_theta,nImages] = size(npf);
nRis_tilde = size(Ris_tilde,3);

% %precompile the shift phases
% shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
% [~,nshifts] = size(shift_phases);

viis = estimate_viis(ciis,Ris_tilde,npf,refq);

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
diffs = zeros(1,nchoosek(nImages,2));
counter = 0;

inds = sub2ind([n_theta/2,n_theta],cijs(:,:,:,:,1),cijs(:,:,:,:,2));
[C,~,IC] = unique(inds(:));

for i=1:nImages
%     if mod(i,10) == 0; reset(g); end;
    npf_i = npf(:,1:n_theta/2,i);
    npf_i(1,:) = 0; % effectivly remove dc term
    % normalize each ray to be of norm 1
    npf_i = bsxfun(@rdivide,npf_i,sqrt(sum((abs(npf_i)).^2)));
    Npf_i = gpuArray(single(npf_i));
    for j=i+1:nImages
        
        t1 = clock;
        
        counter = counter + 1;
        
        npf_j = npf(:,:,j);
        npf_j(1,:) = 0; % effectivly remove dc term
        npf_j = bsxfun(@rdivide,npf_j,sqrt(sum((abs(npf_j)).^2)));
        Npf_j = gpuArray(single(npf_j));
        
        PiPj = Npf_i'*Npf_j;

        Corrs = PiPj(C);
        Corrs = Corrs(IC);
        Corrs = reshape(Corrs,nRis_tilde,nRis_tilde,[],5);
        Corrs = squeeze(mean(real(Corrs),4));
        
        [YY,II] = max(Corrs(:));
        [II,JJ,KK] = ind2sub(size(Corrs),II);
        
        YY = gather(YY);
        II = gather(II);
        JJ = gather(JJ);
        KK = gather(KK);
        
        opt_Rs_tilde(i,j)    = II;
        opt_Rs_tilde(j,i)    = JJ;
        opt_thetaij(i,j)     = KK;
        max_corrs_stats(i,j) = YY;
        
        Ri_tilde = Ris_tilde(:,:,II);
        Rj_tilde = Ris_tilde(:,:,JJ);
        R_theta_ij = R_theta_ijs(:,:,KK);
        
        Rij = Ri_tilde.'*R_theta_ij*Rj_tilde;
        
        Ri_gt = q_to_rot(refq(:,i)).';
        Rj_gt = q_to_rot(refq(:,j)).';
        
        diff = zeros(1,5);
        for s=0:4
            Rij_gs_gt = Ri_gt.'*g^s*Rj_gt;
            diff(s+1) = min(norm(Rij-Rij_gs_gt,'fro'),norm(Rij-J*Rij_gs_gt*J,'fro'));
        end
        diffs(counter) = min(diff);
        
        vijs(:,:,counter) = 1/5*Ri_tilde.'*(R_theta_ij+g*R_theta_ij+g^2*R_theta_ij+g^3*R_theta_ij+g^4*R_theta_ij)*Rj_tilde;
        
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

mse_vij = sum(diffs.^2)/numel(diffs);
log_message('MSE of vij: %e',mse_vij);

end