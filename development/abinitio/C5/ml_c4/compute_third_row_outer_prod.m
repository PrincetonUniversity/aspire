function [vijs,viis,max_corrs_stats] = ...
    compute_third_row_outer_prod(npf,ciis,cijs,Ris_tilde,R_theta_ijs,max_shift,shift_step,is_handle_equators,refq)

[n_r,n_theta,nImages] = size(npf);
nRis_tilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

%precompile the shift phases
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
viis = estimate_viis(ciis,Ris_tilde,npf,shift_phases,is_handle_equators,refq);

g_shift_phases = gpuArray(double(shift_phases));
[~,nshifts] = size(shift_phases);

opt_Rs_tilde    = zeros(nImages,nImages);
opt_thetaij     = zeros(nImages,nImages);
max_corrs_stats = zeros(nImages,nImages);

log_message('\ncomputing likelihood');
msg = [];

vijs = zeros(3,3,nchoosek(nImages,2));
g = [0 -1 0; 1 0 0; 0 0 1]; % a rotation of 90 degrees about the z-axis
J = diag([1,1,-1]);
% g = gpuDevice(1);
diffs = zeros(1,nchoosek(nImages,2));
counter = 0;

inds = sub2ind([n_theta,n_theta/2],cijs(:,:,:,:,1),cijs(:,:,:,:,2));
[C,~,IC] = unique(inds(:));

for i=1:nImages
    %     if mod(i,10) == 0; reset(g); end;
    npf_i = npf(:,:,i);
    g_npf_i = gpuArray(double(npf_i));
    % ignoring dc-term
    g_npf_i(1,:) = 0;
    
    % nomalize each ray to be norm 1
    norms   = sqrt(sum((abs(g_npf_i)).^2));
    g_npf_i = bsxfun(@rdivide,g_npf_i,norms);
    
    for j=i+1:nImages
        
        t1 = clock;
        
        counter = counter + 1;
        
        npf_j = npf(:,1:n_theta/2,j);
        g_npf_j = gpuArray(double(npf_j));
        
        Corrs = zeros([numel(C),1],'gpuArray');
        for s=1:nshifts
            g_npf_j_shifted = bsxfun(@times,g_npf_j,g_shift_phases(:,s));
            
            % ignoring dc-term
            g_npf_j_shifted(1,:) = 0;
            
            % nomalize each ray to be norm 1
            norms           = sqrt(sum((abs(g_npf_j_shifted)).^2));
            g_npf_j_shifted = bsxfun(@rdivide,g_npf_j_shifted,norms);
            
            
            PiPj = g_npf_i'*g_npf_j_shifted;
            
            Corrs_s = PiPj(C);
            
            Corrs = max([Corrs real(Corrs_s(:))],[],2);
        end
        
        Corrs = Corrs(IC);
        Corrs = reshape(Corrs,[nRis_tilde,nRis_tilde,n_theta_ij/4,4]);
        Corrs = squeeze(mean(real(Corrs),4));
        
        [YY,II] = max(Corrs(:));
        
        [II,JJ,KK] = ind2sub([nRis_tilde,nRis_tilde,n_theta_ij/4],II);
        
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
        
        % TODO: move outside the loop
        if exist('refq','var') && ~isempty(refq)
            
            
            Rij = Ri_tilde.'*R_theta_ij*Rj_tilde;
            
            Ri_gt = q_to_rot(refq(:,i)).';
            Rj_gt = q_to_rot(refq(:,j)).';
            
            diff = zeros(1,4);
            for s=0:3
                Rij_gs_gt = Ri_gt.'*g^s*Rj_gt;
                diff(s+1) = min(norm(Rij-Rij_gs_gt,'fro'),norm(Rij-J*Rij_gs_gt*J,'fro'));
            end
            diffs(counter) = min(diff);
        end
        
        vijs(:,:,counter) = 0.25*Ri_tilde.'*(R_theta_ij+g*R_theta_ij+g^2*R_theta_ij+g^3*R_theta_ij)*Rj_tilde;
        
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
    mse_vij = sum(diffs.^2)/numel(diffs);
    log_message('MSE of vij: %e',mse_vij);
end

end