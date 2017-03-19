function [rots,in_plane_rotations] = estimate_inplane_rotations_nogpu(npf,vis,inplane_rot_res,max_shift,shift_step)

% if params.isRemoveNonRank1
%     W = conf;
% elseif params.isUseWeight_in_plain_rot
%     W = conf + conf' ; %+ eye(K);
%     D = sum(W, 2);
%     W = diag(D)\W;  % normalize every row
%     % adding 'I' for debug purposes. Namely that in clean setting thw spectrum is indeed (1,0,0,...0).
%     %Otherwise everything is shifted by 1
%     W = W + eye(nImages);
% else
%     W = ones(nImages,nImages);
% end

[n_r,n_theta,nImages] = size(npf);

shift_phases = calc_shift_phases(n_r,max_shift,shift_step);

[~,nshifts] = size(shift_phases);
assert(nImages == size(Ris_tilde,3));

H = zeros(nImages,nImages);

%TODO: see if params.inplane_rot_res=0.5 gives better results than, say, 1
%(and if so consider also 0.25 resolution)
resolution = inplane_rot_res; % assumed to be in degrees
theta_ij = (0:resolution:(360-resolution))*pi/180;
n_theta_ij = numel(theta_ij);

assert(mod(n_theta_ij,4)==0);

cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ij = [cos_theta_ij  ; sin_theta_ij; zrs; ...
    -sin_theta_ij ; cos_theta_ij; zrs; ...
    zrs;            zrs;          ons];

R_theta_ij = reshape(R_theta_ij,3,3,n_theta_ij);

max_corrs = zeros(1,nchoosek(nImages,2));
max_idx_corrs = zeros(1,nchoosek(nImages,2));
counter = 0;

% npf = filter_npf(npf);
%TODO : remove normalizition as this should be done after we get the
%shifted image
% norms = sqrt(sum((abs(npf)).^2));
% npf = bsxfun(@rdivide,npf,norms);

Ris_tilde = zeros(3,3,nImages);
for i = 1:nImages
    v_i = vis(:,i).';
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(v_i);
end

log_message('ignoring dc-term');
npf(1,:,:) = 0;

npf_normalized = normalize_ray(npf);

npairs = nchoosek(nImages,2);
fprintf('\n');printProgressBarHeader;
for i = 1:nImages
    
    npf_i = npf(:,:,i);
    
    npf_i_shifted = get_shifted_cps(npf_i,shift_phases);
    %     npf_i_shifted = reshape(npf_i_shifted,n_r,n_theta*nshifts);
    
    npf_i_shifted = normalize_ray(npf_i_shifted);
%     Npf_i_shifted = gpuArray(double(npf_i_shifted));
    
    for j = i+1:nImages
        
        counter = counter+1;
        
        progressTic(counter,npairs);
        
        Ri_tilde = Ris_tilde(:,:,i);
        Rj_tilde = Ris_tilde(:,:,j);
        
        
        npf_j = npf_normalized(:,:,j);
        
        tmp = multiprod(Ri_tilde.',R_theta_ij);
        Us  = multiprod(tmp,Rj_tilde);
        
        Us = reshape(Us,9,n_theta_ij);
        % extract a common-line index for each possible theta_ij
        c_1 = [-Us(8,:) ;  Us(7,:)]; %[-U(2,3)  U(1,3)];
        c_2 = [ Us(6,:) ; -Us(3,:)]; %[ U(3,2) -U(3,1)];
        
        cijs = clAngles2Ind(c_1,n_theta);
        cjis = clAngles2Ind(c_2,n_theta);
        
%         Npf_j = gpuArray(double(npf_j));
        
        co = bsxfun(@times,npf_i_shifted(:,cijs,:),conj(npf_j(:,cjis))); % cross-correaltion - so we want to conjugate
        corrs = sum(co);
%         corrs = gather(Corrs);
        %         corrs = sum(conj(npf_i(:,cijs)).*npf_j(:,cjis)); % cross-correaltion - so we want to conjugate
        
        corrs = reshape(corrs,n_theta_ij/4,4,nshifts);
        
        if nshifts > 1
            % the 1d shifts depend not only on the 2d shifts of each and every
            % image, but also on their common-lines. That is, for a given image, the 1d
            % shift might be different for different lines !
            % In particular, this means that each one of the 4 lines we are considering in a given image
            % may correspond to a DIFFERENT 1d shift. Therefore each line has the dof, so to speak, of choosing its own shift.
            % All in all, this means that we max out the shift BEFORE taking
            % the mean score over all 4 lines.
            corrs = max(real(corrs),[],3);
        end
        % now take the mean score over all 4 lines, and find the maximum
        % quadraple of lines
        [max_corr,max_idx_corr] = max(mean(real(corrs),2));
        max_corrs(counter) = max_corr; % this is only for stats
        max_idx_corrs(counter) = max_idx_corr; % this is only for stats
        theta_diff = resolution*(max_idx_corr-1)*pi/180;
        
        H(i,j) = cos(4*theta_diff) + sqrt(-1)*sin(4*theta_diff);
    end
end
% note entry (i,j) corresponds to exp^(-i(-theta_i+theta_i)). Therefore to
% construct the hermitian matrix, entry (j,i) is the **conjugate** of entry
% (i,j)
% transpose
H = H + H' + eye(nImages); % put 1 on diagonal since : exp^(-i*0) = 1

[v, d] = eigs(H, 10, 'lm');
[evals, ind] = sort(diag(d), 'descend');
evect1 = v(:,ind(1)); %zi = exp(-i*2*ai), Ri = Rot(ai)*Ri0

disp(' ');
disp(['First 5 eigenvalues : ' num2str(evals(1:5)')]);

% figure,
% hist(evals, 40);
% set(gca, 'FontSize', 20);
% title('SO(2) Sync');

in_plane_rotations = zeros(3,3,nImages);
for i  = 1:nImages
    zi  = evect1(i);
    zi  = zi/abs(zi); % rescale so it lies on unit circle
    c   = real(zi^(1/4));  % we have four time the desired angle
    s   = -imag(zi^(1/4)); % we have four time the desired angle
    in_plane_rotations(:,:,i) = [c -s 0; s c 0; 0 0 1];
end

rots = form_rotation_matrices(Ris_tilde,in_plane_rotations);


end


function R = complete_3rdRow_to_rot(r3)

% build R \in SO(3) from the 3rd Row
% 3rd col is not aligning to z-axis

tmp = sqrt(r3(1)^2 + r3(2)^2);

r1 = [r3(2)/tmp -r3(1)/tmp 0];
r2 = [r3(1)*r3(3)/tmp r3(2)*r3(3)/tmp -tmp];

R = [r1; r2; r3];

end




% function in_plane_rotations = test_estimate_inplane_rotations(npf,shift_phases,params)
% 
% nImages = params.K;
% Ris_tilde_gt = zeros(3,3,nImages);
% 
% for i=1:nImages
%     
%     Ri = q_to_rot(params.refq(:,i))';
%     qi_gt = Ri(3,:);
%     Ris_tilde_gt(:,:,i) = complete_3rdRow_to_rot(qi_gt);
% end
% 
% in_plane_rotations = estimate_inplane_rotations(npf,Ris_tilde_gt,shift_phases,params);
% 
% 
% end