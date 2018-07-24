function all_self_corrs = estimate_viis(ciis,Ris_tilde,npf,n_symm,shift_phases,refq)

% if is_handle_equators && mod(n_symm,2) ~= 0
%     error('can only handle equators if symmetry class is even');
% end

[~,n_theta,nImages] = size(npf);

% find which of the candidates rotations are equator images
cii_equators_inds = squeeze(find(abs(acosd(Ris_tilde(3,3,:)) - 90) < 10));
% if is_handle_equators && ~isempty(cii_equators_inds)
%     
%     if ~exist('refq','var'), refq = []; end;
%     
%     [viis_equators,inds_eq_images] = ...
%         handle_equator_images(npf,Ris_tilde,ciis,cii_equators_inds,max_shift,shift_step,10,0.1,refq);
%     
%     viis(:,:,inds_eq_images) = viis_equators;
%     
% else
%     inds_eq_images = [];
% end

g_shift_phases = gpuArray(double(shift_phases));
[~,nshifts] = size(shift_phases);

% there are alltogether n_symm-1 self-cl in each image. If n_symm is even then the pairs of one of these
% lines (the one in the middle) is itself and is therefore useless (it is the one that lies on the global X_Y plane)
[n_selfcl_pairs_times2, ncands] = size(ciis);

n_selfcl_pairs = n_selfcl_pairs_times2/2;
assert(n_selfcl_pairs == floor((n_symm-1)/2)); 

inds_scls = zeros(ncands,n_selfcl_pairs); 

for p=1:n_selfcl_pairs
    inds_scls(:,p) = sub2ind([n_theta,n_theta/2],ciis(2*p-1,:),ciis(2*p,:)); 
end

opt_cand_ind_statis = zeros(1,nImages); % just for statistics the optimal candidate per image
max_corrs_stats     = zeros(1,nImages); % just for statistics

all_self_corrs = zeros(ncands,nImages); % for each image gives the correlation of all pairs of self common-lines induced by each candidate

for i=1:nImages
    % nomalize each ray to be norm 1
    npf_i = npf(:,:,i);
    norms   = sqrt(sum((abs(npf_i)).^2));
    npf_i = bsxfun(@rdivide,npf_i,norms);
    npf(:,:,i) = npf_i;
end

msg = [];
for i=1:nImages
    
    t1 = clock;
%     if ismember(i,inds_eq_images); 
%         continue; 
%     end
    npf_i = npf(:,:,i);
    g_npf_i = gpuArray(double(npf_i));
    
    
    g_half_npf_i = g_npf_i(:,1:n_theta/2);
    
    Corrs_scls = zeros([ncands,n_selfcl_pairs],'gpuArray');

    for s=1:nshifts
        g_half_npf_i_shifted = bsxfun(@times,g_half_npf_i,g_shift_phases(:,s));
        
        
        PiPi = g_npf_i'*g_half_npf_i_shifted; % compute the correlation of Pi with itself for a given shift
        
        for p=1:n_selfcl_pairs
            
            Corrs_scls_tmp = PiPi(inds_scls(:,p));
            % the optimal shift in fourier space is per line and not per
            % image. Thus let each line "choose" its optimal shift
            Corrs_scls(:,p) = max([Corrs_scls(:,p) real(Corrs_scls_tmp(:))],[],2);
        end
    end
    
    
    Corrs = mean(Corrs_scls,2);
    
    % equaotr candidates induce collinear self common-lines which always have a perfect correlation. Hence we must exclude these
    Corrs(cii_equators_inds) = -inf; 
    
    Corrs = real(Corrs);
    [Y,I] = max(Corrs);
   
    % just for statistics
    opt_cand_ind_statis(i) = gather(I);
    max_corrs_stats(i)   = gather(Y);
    
    % output variable
    all_self_corrs(:,i) = gather(Corrs);
    
    if true
        t2 = clock;
        t = etime(t2,t1);
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('i=%3d/%3d t=%7.5f',i,nImages,t);
        fprintf('%s',msg);
    end
        
end

% be polite and return correlation zero instead of infinity for all equator
% candidates
all_self_corrs(cii_equators_inds,:) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%!!!!!!!!!!!! DEBUG CODE !!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('refq','var') && ~isempty(refq)
    debug_analyze_viis(opt_cand_ind_statis, n_symm, Ris_tilde, n_theta, refq)
end

end

function debug_analyze_viis(opt_cand_ind_statis, n_symm, Ris_tilde, n_theta, refq)

g = [cosd(360/n_symm) -sind(360/n_symm) 0;
    sind(360/n_symm)  cosd(360/n_symm) 0;
    0 				 0  1]; % a rotation of 360/n_symm degrees about the z-axis
J = diag([1,1,-1]);

if exist('refq','var') && ~isempty(refq)
    
    nImages = size(refq,2);
    viis = zeros(3,3,nImages);
    for i=1:nImages
        
        Ri_opt = Ris_tilde(:,:,opt_cand_ind_statis(i));
        
        sum_Rii = zeros(3,3);
        for s=0:n_symm-1
            sum_Rii = sum_Rii + Ri_opt.'*g^s*Ri_opt;
        end
        viis(:,:,i) = sum_Rii/n_symm;
        
    end
    %% check mse of vij
    errs_mse = zeros(1,nImages);
    for i=1:nImages
        
        rot_i_gt = q_to_rot(refq(:,i)).';
        
        vii_gt = rot_i_gt(3,:).'* rot_i_gt(3,:);
        
        vii = viis(:,:,i);
        
        err = min(norm(vii-vii_gt,'fro'),norm(J*vii*J-vii_gt,'fro'));
        errs_mse(i) = err;
        
    end
    
    log_message('\n vii errs mse=%.2f',mean(errs_mse));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    diffs = zeros(1,nImages);
    for i=1:nImages
        Ri_gt = q_to_rot(refq(:,i)).';
        vi_gt = Ri_gt(3,:);
        vii_gt = vi_gt.'*vi_gt;
        
        vii = viis(:,:,i);
        
        diffs(i) =  min([norm(  vii     - vii_gt,'fro'),...
            norm(  vii.'   - vii_gt,'fro')...
            norm(J*vii*J   - vii_gt,'fro'),...
            norm(J*vii.'*J - vii_gt,'fro')]);
    end
    
    mse_vii = sum(diffs.^2)/numel(diffs);
    log_message('MSE of vii: %e',mse_vii);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    is_correct = zeros(1,nImages);
    hand_idx = zeros(1,nImages);
    angle_tol_err = 10/180*pi; % how much angular deviation we allow for a self-common-line to have
    for i=1:nImages
            
        Ri_gt = q_to_rot(refq(:,i)).';
        
        
        Rii_g_gt = Ri_gt.'*g*Ri_gt;
        
        c_1_gt = [-Rii_g_gt(8) ;  Rii_g_gt(7)]; %[-Riig(2,3)  Riig(1,3)];
        c_2_gt = [ Rii_g_gt(6) ; -Rii_g_gt(3)]; %[ Riig(3,2) -Riig(3,1)];
        
        cigi_gt = clAngles2Ind(c_1_gt,n_theta);
        cgii_gt = clAngles2Ind(c_2_gt,n_theta);
        
        
        
        Ri_tilde   = Ris_tilde(:,:,opt_cand_ind_statis(i));
        
        diff_s = zeros(1,n_symm-1);
        hand_idx_s = zeros(1,n_symm-1);
        for s=1:11
            Rii_g = Ri_tilde.'*g^(s)*Ri_tilde;
            
            c_1 = [-Rii_g(8) ;  Rii_g(7)]; %[-Riig(2,3)  Riig(1,3)];
            c_2 = [ Rii_g(6) ; -Rii_g(3)]; %[ Riig(3,2) -Riig(3,1)];
            
            cigi = clAngles2Ind(c_1,n_theta);
            cgii = clAngles2Ind(c_2,n_theta);
            
            
            cigi_diff = (cigi_gt-cigi)*2*pi./n_theta;
            cgii_diff = (cgii_gt-cgii)*2*pi./n_theta;
            
            % take absolute cosine because of handedness
            % there might be +180 independendt diff for each image which at this stage
            % hasn't been taken care yet.
            diff = acos(cos(cigi_diff))    + acos(cos(cgii_diff));
            diff_J = acos(cos(cigi_diff+pi)) + acos(cos(cgii_diff+pi));
            if diff < diff_J
                diff_s(s) = diff;
                hand_idx_s(s) = 1;
            else
                diff_s(s) = diff_J;
                hand_idx_s(s) = 2;
            end
        end
        [min_diff_s,min_idx] = min(diff_s);
        if min_diff_s < 2*angle_tol_err
            is_correct(i) = 1;
            hand_idx(i) = hand_idx_s(min_idx);
        end
    end
    
    scl_dist = histc(hand_idx,1:2)/nImages;
    scl_detec_rate = sum(is_correct)/nImages;
    log_message('\nself common lines detection rate=%.2f%%',scl_detec_rate*100);
    log_message('scl_J_dist=[%.2f %.2f]',scl_dist);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bad_inds = find(diffs > 0.1);
    bad_polar_angs = zeros(1,numel(bad_inds));
    for i=1:numel(bad_inds)
        ind = bad_inds(i);
        bad_rot = q_to_rot(refq(:,ind)).';
        bad_polar_angs(i) = acosd(bad_rot(3,3));
    end
    figure; hist(bad_polar_angs,180);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%!!!!!!!!!!!! END OF DEBUG CODE !!!!!!!!!!!!!!%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



function inds_eq_images = detect_equator_images(npf,max_shift,res_factor,fraction,refq)
%
% Finds the images who's corresponding beaming direction is close to
% (*,#,eps)^T where eps is a small number and *,# are any two numbers. It
% is based on the fact that equator images of a C4 symmetric molecule
% have a reflection symmetry about the horizontal axis. That is, im(:,theta) = im(:,-theta)
%
%
% Input parameters:
%   npf             A 3d table containing the 2-d fourier transform of each
%                   projection image (nImages is the size of the third dimension)
%   max_shift       (Optional) The maximum spatial shift that each image is
%                   assumed to have. Default:15
%   res_factor      (Optional) Angular resolution factor that each image
%                   should undergo. For example if res_factor=10 then all
%                   values along angles 0-9 are concatenated (and same for 10-19,20-29, etc)
%   fraction        (Optional) the fraction of input images
%                   to decalare as equator images. Defualt=0.1
%   refq            (Optional) A 2d table where the i-th column is the
%                   quaternion that corresponds to the beaming direction of
%                   the i-th image.
%
%
% Output parameters:
%   inds_eq_images  The indexes of equator images found. The number of
%                   indexes is fraction*nImages
%

log_message('Detecting equator images');

[n_r,n_theta,nImages] = size(npf);

nshifts = (2*max_shift+1)^2;
shifted_deltas = zeros(n_r,n_r,nshifts);
i_shift = 1;
for s_x = -max_shift:max_shift
    for s_y  = -max_shift:max_shift
        shifted_deltas(ceil(n_r/2)+s_x,...
            ceil(n_r/2)+s_y,...
            i_shift) = 1;
        i_shift = i_shift + 1;
    end
end
[phases,~] = cryo_pft(shifted_deltas, n_r, n_theta, 'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 3: detect the equotor images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

images_score = zeros(1,nImages);
Phases = gpuArray(single(phases));
msg = [];

for i_img = 1:nImages
    
    t1 = clock;
    
    pim = npf(:,:,i_img);
    g_pim = gpuArray(single(pim));
    g_pim_shifted = bsxfun(@times,g_pim,Phases);
    
    g_norms = sqrt(sum((abs(g_pim_shifted)).^2,2));
    g_pim_shifted = bsxfun(@rdivide,g_pim_shifted,g_norms);
    %     pim_shifted = gather(g_pim_shifted);
    
    g_pim_shifted(1,:,:) = 0; %ignore the dc term;
    
    % Equator images of a C4 symmetric molecule have the proprty that each
    % of its images has a reflection symmetry about the horizontal axis. That is,
    % im(:,theta) = im(:,-theta)
    
    %flip all images with respect to theta
    g_pim_shifted_flipped = flipdim(g_pim_shifted,2);
    
    %     % to find the translation we compute the cross power spectrum
    %     g_pim_shifted            = gpuArray(single(pim_shifted));
    %     Pim_shifted_flipped      = gpuArray(single(pim_shifted_flipped));
    
    Cross_pwr_spec = fft(g_pim_shifted ,[],2).*conj(fft(g_pim_shifted_flipped,[],2));
    Inv_cross_pwr_spec = ifft(Cross_pwr_spec,[],2);
    inv_cross_pwr_spec = gather(Inv_cross_pwr_spec);
    
    [nr,nc,~] = size(inv_cross_pwr_spec);
    inv_cross_pwr_spec = reshape(inv_cross_pwr_spec,nr*res_factor,nc/res_factor,[]);
    [shifts_score,~] = max(real(sum(inv_cross_pwr_spec,1)));
    images_score(i_img) = max(shifts_score);
    
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t2 = clock;
    t = etime(t2,t1);
    bs = char(repmat(8,1,numel(msg)));
    fprintf('%s',bs);
    msg = sprintf('k=%3d/%3d t=%7.5f',i_img,nImages,t);
    fprintf('%s',msg);
    
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

[~, sorted_inds] = sort(images_score,'descend');

inds_eq_images = sorted_inds(1:floor(fraction*numel(sorted_inds)));


%%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check how well we did in detecting equator images
if exist('refq','var') && ~isempty(refq)
    assert(nImages == size(refq,2));
    is_eq_gt = false(1,nImages);
    for k = 1:nImages
        Rk_gt  = q_to_rot(refq(:,k))';
        if( abs(Rk_gt(3,3)) < cosd(85))
            is_eq_gt(k) = true;
        end
    end
    TPR = 100*sum(is_eq_gt(inds_eq_images))/sum(is_eq_gt);
    log_message('True-positive-rate of equator images=%.2f%% (#detected_equators/#total_equators)',TPR)
    %     figure; plot(vv(is_eq),'g*'); hold on; plot(vv(~is_eq),'r*');
end

%%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if ismember(nImages,inds_eq_images)
%     printf('\nTop-view image was accidently detected as equator image.');
%     inds_eq_images = inds_eq_images(inds_eq_images~=params.K);
% end

% printf('\nRemoving %d images (%.2f%%) in order eliminate equator images',numel(inds_eq_images),numel(inds_eq_images)/params.K*100);
% masked_projs(:,:,inds_eq_images) = [];
% noisy_projs(:,:,inds_eq_images) = [];
% npf(:,:,inds_eq_images) = [];
% params.K = params.K - numel(find(inds_eq_images));
% if ~params.real_data
%     params.refq(:,inds_eq_images) = [];
% end

end


function [viis_equators,inds_eq_images] = ...
    handle_equator_images(npf,Ris_tilde,ciis,cii_equators_inds,max_shift,shift_step,res_factor,fraction,refq)

if ~exist('fraction','var')
    fraction = 0.1;
end

if ~exist('res_factor','var')
    res_factor = 10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 1: detect equator images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('refq','var'), refq = []; end;
inds_eq_images = detect_equator_images(npf,max_shift,res_factor,fraction,refq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 2: detect self-common-lines for all equator images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_message('Detecting self-common-lines for equator images');
% Problem: the angles between the first and third self-cl for equaot images
% is 180 degrees. As these are conjugate symmetric for any image, these
% cannot be detected.
% Observation1: ALL equator images intersect at the global z-axis. So do in
% particualr any pair of equator images Ri and gRi.
% Observation2: each pair of equator images have two pairs of common-lines
% (not four)
% Solution: find the self-common-lines for each equator image by searching the
% median common-line with all other equator images
%
% clmatrix_eq = commonlines_gaussian(npf_eq,params.max_shift,1); %TODO: extract shift info for debug
nEquator_images = numel(inds_eq_images);
%     scl_equators = zeros(2,nEquators);
clmatrix_eq = cryo_clmatrix_gpu(npf(:,:,inds_eq_images),nEquator_images,0,max_shift,shift_step);
% we'll take the median line, so make sure all lines are in [0,180]
n_theta = size(npf,2);
clmatrix_eq = mod(clmatrix_eq-1,n_theta/2)+1;
scl_equators = zeros(2,nEquator_images);
scl_equators(1,:) = median(clmatrix_eq,2)';
scl_equators(2,:) = mod(scl_equators(1,:)+n_theta/2-1,n_theta)+1; % we **know** they are 180 degrees apart.

ciis_equators = ciis(:,cii_equators_inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 4: assign to each equator image an equator rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
viis_equators = zeros(3,3,nEquator_images);
g = [cosd(360/n_symm) -sind(360/n_symm) 0; 
	 sind(360/n_symm)  cosd(360/n_symm) 0; 
	 0 				 0  1]; % a rotation of 360/12 degrees about the z-axis
for i=1:nEquator_images
    % it is fine to get the self-common-lines in the wrong order
    [YYa,IIa] = min(sum(abs(bsxfun(@minus,scl_equators(:,i),       double(ciis_equators)))));
    [YYb,IIb] = min(sum(abs(bsxfun(@minus,scl_equators(:,i),flipud(double(ciis_equators))))));
    if YYa < YYb
        Ri_opt = Ris_tilde(:,:,IIa);
    else
        Ri_opt = Ris_tilde(:,:,IIb);
    end
    Rii = Ri_opt.'*g*Ri_opt;
    viis_equators(:,:,i) = 0.5*(Rii+Rii.');
end

end