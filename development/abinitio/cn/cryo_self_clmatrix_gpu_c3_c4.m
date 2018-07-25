function [sclmatrix,correlations,shifts] = cryo_self_clmatrix_gpu_c3_c4(n_symm,npf,max_shift,shift_step,is_handle_equator_ims,refq)
% Input parameters:

%   n_symm                  Either 3 (for c_3) or 4 (for c_4)
%   npf                     A 3D array where each image npf(:,:,i) corresponds to the Fourier
%                           transform of projection i.
%   max_shift               The maximum spatial shift that each image is
%                           assumed to have. Default:15
%   shift_step              (Optional) Default:0.5
%   is_handle_equator_ims   (Optional) whether to handle equator images 
%                           seperatelly or not. Defualt:true 
%   equator_res_fact        (Optional) Angular resolution factor that each image
%                           should undergo. For example if res_factor=10 then all
%                           values along angles 0-9 are concatenated (and same for 10-19,20-29, etc)
%   equator_fraction        (Optional) the fraction of input images to decalare as
%                           equator images. Default:0.1
%   refq                    (Optional) A 2d table where the i-th column is the
%                           quaternion that corresponds to the beaming direction of
%                           the i-th image.
% Output parameters:
%   sclmatrix               A 2xnImages table where the i-th column holds
%                           the indexes of the first and third
%                           self-common-line in image i. 
%   correlations            An array of length nImages whose i-th entry
%                           holds the correlation between the two self common-lines that were found
%                           in image i
%   shifts                  An array of length nImages whose i-th entry
%                           holds the shift found for image i

log_message('detecting self-common-lines');

if n_symm ~= 3 && n_symm ~= 4
    error('n_symm may be either 3 or 4');
end

if ~exist('equator_fraction','var')
    equator_fraction = 0.1;
end

if ~exist('equator_res_fact','var')
    equator_res_fact = 10;
end

if ~exist('is_handle_equator_ims','var')
    if n_symm == 4
        is_handle_equator_ims = true;
    else
        is_handle_equator_ims = false;
    end
end

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 15;
end

log_message('detecting self-common-lines');

[n_r,n_theta,nImages] = size(npf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: detect self-common-lines for all images
% Step 2: detect all equator images (only for c4)
% Step 3: detect self-common-lines for all equator images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for the angle between self-common-lines. In theory it is
% [60,180] but all 180 apart lines are perfectly correlated so we set it to
% be smaller.
if n_symm == 3
    min_angle_diff = 60*pi/180;
    max_angle_diff = 165*pi/180;
else % i.e. n_symm == 4
    min_angle_diff = 90*pi/180;
    max_angle_diff = 160*pi/180;
end

% the self-common-line matrix holds per image two indeces that represent
% the two self common-lines in the image
sclmatrix = zeros(2,nImages);
correlations = zeros(1,nImages);
shifts = zeros(1,nImages);

% the angular difference between each two self common lines in a given
% image is [90,180], so create a mask.
[X,Y] = meshgrid(1:n_theta/2,1:n_theta);
diff = Y-X;
unsigned_angle_diff = acos(cos(diff.*2*pi./n_theta));

good_diffs = unsigned_angle_diff > min_angle_diff & ...
    unsigned_angle_diff < max_angle_diff;

shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
[n_r_,nshifts] = size(shift_phases);
assert(n_r == n_r_);
g_shift_phases = gpuArray(single(shift_phases));

msg = [];
for i=1:nImages
    
    t1 = clock;
    
    npf_i = npf(:,:,i);
    % each image is conjugate symmetric so need only consider half of the lines
    pi_half = npf_i(:,1:n_theta/2);
    
    g_npf_i   = gpuArray(single(npf_i));
    g_pi_half = gpuArray(single(pi_half));
    
    % generate all shifted copies of the image
    g_pi_half_shifted = zeros([size(pi_half),nshifts],'gpuArray');
    for s=1:nshifts
        g_pi_half_shifted(:,:,s) = bsxfun(@times,g_pi_half,g_shift_phases(:,s));
    end
    
    g_pi_half_shifted = reshape(g_pi_half_shifted,n_r,n_theta/2*nshifts);
    
    % ignoring dc-term
    g_npf_i(1,:) = 0;
    g_pi_half_shifted(1,:) = 0;
    
    % nomalize each ray to be norm 1
    norms = sqrt(sum((abs(g_npf_i)).^2));
    g_npf_i = bsxfun(@rdivide,g_npf_i,norms);
    
    % nomalize each ray to be norm 1
    norms = sqrt(sum((abs(g_pi_half_shifted)).^2));
    g_pi_half_shifted = bsxfun(@rdivide,g_pi_half_shifted,norms);
    
    Corr = g_npf_i.'*g_pi_half_shifted;
    corr = gather(Corr);
    corr = reshape(corr,[n_theta,n_theta/2,nshifts]);
    
    corr = bsxfun(@times,corr,good_diffs);
    
    [correlation,idx] = max(real(corr(:)));
    
    [sc1, scl3, shift] = ind2sub([n_theta, n_theta/2, nshifts],idx);
    sclmatrix(:,i) = [sc1, scl3]';
    correlations(i) = correlation;
    shifts(i) = shift; %TODO: need to "translate" shift 2*max_shift+1 = ...
    
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    t = etime(t2,t1);
    bs = char(repmat(8,1,numel(msg)));
    fprintf('%s',bs);
    msg = sprintf('k=%3d/%3d t=%7.5f',i,nImages,t);
    fprintf('%s',msg);
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

log_message('\nself-common-lines correlation (median)=%.2f',median(correlations));

if is_handle_equator_ims
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 2  : detect equator images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inds_eq_images = detect_equator_images(npf,max_shift,equator_res_fact,...
        equator_fraction);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Step 3: detect self-common-lines for all equator images
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
    nEquators = numel(inds_eq_images);
%     scl_equators = zeros(2,nEquators);
    clmatrix_eq = cryo_clmatrix_gpu(npf(:,:,inds_eq_images),...
                                    nEquators,0,max_shift,shift_step);
    % we'll take the median line, so make sure all lines are in [0,180]
    clmatrix_eq = mod(clmatrix_eq-1,n_theta/2)+1; 
    sclmatrix(1,inds_eq_images) = median(clmatrix_eq,2)';
    sclmatrix(2,inds_eq_images) = ...
        mod(sclmatrix(1,inds_eq_images)+n_theta/2-1,n_theta)+1; % we **know** they are 180 degrees apart.
end

% analysis against ground-truth
if exist('refq','var') && ~isempty(refq)
    scl_detection_rate(n_symm,sclmatrix,n_theta,refq);
end

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

if ~exist('removal_frac','var')
    fraction = 0.1;
end

if ~exist('res_factor','var')
    res_factor = 10;
end

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




function detec_rate = scl_detection_rate(n_symm,sclmatrix,n_theta,refq)
%
% Calculates the detection rate of self-common-lines. It is invariant to a
% possible J-ambiguity (independently for each image), as well as invariant
% to an RigRi<->Rig^{3}Ri ambiguity (independently for each image)
% 
% Input parameters:
%   sclmatrix  A 2xnImages table where the i-th column holds
%              the indexes of the first and third
%              self-common-line in image i. 
%   n_theta    The angular resolution of the images. E.g.,n_theta=360 
%              means that there are 360 lines in each images 
%   refq       A 4-by-n table. The i-th column represent the quaternion of
%              that corresponds to the rotation matrix of the i-th image
%
% Output parameters:
%   detec_rate The detection rate in [0,1] of self common-lines against the
%               ground truth
              
angle_tol_err = 10/180*pi;
% Two issues:
% 1. DOF: cannot tell the difference between first\third self-common-line
% 2. Ambiguity: handedness - At this stage the handedness is independent for each image
nImages = size(sclmatrix,2);
% sclmatrix_correct = zeros(size(sclmatrix));
sclmatrix_gt = detectScls_gt(n_symm,n_theta,refq); % clmatrix_gt is a 2*n matrix

sclmatrix_diff1 = sclmatrix_gt - sclmatrix;
sclmatrix_diff2 = sclmatrix_gt - flipud(sclmatrix); % we cannot (and need not) tell the difference between the first scl and third scl
clmatrix_diff1_angle = sclmatrix_diff1*2*pi./n_theta;
clmatrix_diff2_angle = sclmatrix_diff2*2*pi./n_theta;

% 1. cosine is invariant to 2pi.
% 2. abs is invariant to +-pi diff corresponding to J-ambiguity
clmatrix_diff1_angle_mean = mean(acos(abs(cos(clmatrix_diff1_angle))));
clmatrix_diff2_angle_mean = mean(acos(abs(cos(clmatrix_diff2_angle))));

% we need not tell the difference between the first scl and third scl so
% take the min difference
[min_mean_angle_diff, scl_idx] = min([clmatrix_diff1_angle_mean ; clmatrix_diff2_angle_mean]);
correct_idxs = find(min_mean_angle_diff < angle_tol_err);
detec_rate = numel(correct_idxs)/nImages;
% for debug purposes: just making sure that more or less half of time we
% retrieved the first self-cl and half of the time the third self-cl
scl_dist = histc(scl_idx,1:2)/numel(scl_idx);
log_message('self-common lines detection rate=%.2f%%',detec_rate*100);
log_message('scl-distribution=[%.2f %.2f]',scl_dist);


% find the polar angles of viewing directions of offending rotations
rots_gt = zeros(3,3,nImages);
for k = 1:nImages
    rots_gt(:,:,k) = q_to_rot(refq(:,k))';
end

bad_idxs  = find(min_mean_angle_diff >= angle_tol_err);
if ~isempty(bad_idxs)
    bad_polar_angs = acosd(rots_gt(3,3,bad_idxs));
    nbad = min(numel(bad_polar_angs),5);
    log_message('sample of polar angles failing self-cl [%.2f,%.2f,%.2f,%.2f,%.2f]',...
        bad_polar_angs(1:nbad));
end

end

function sclmatrix_gt = detectScls_gt(n_symm,n_theta,refq)

if n_symm ~= 3 && n_symm ~= 4
    error('n_symm may be either 3 or 4');
end

nImages = size(refq,2);
% we find the first and last self common-lines (c3 and c4s ymmetry)
sclmatrix_gt = zeros(2,nImages);
% correlations_selfcl  = zeros(1,nImages);
rots_gt = zeros(3,3,nImages);
for i = 1:nImages
    rots_gt(:,:,i) = q_to_rot(refq(:,i))';
end

g = [cosd(360/n_symm) -sind(360/n_symm) 0; ...
     sind(360/n_symm)  cosd(360/n_symm) 0; ...
     0                 0  1]; % rotation matrix of 120 or 90 degress around z-axis


for i=1:nImages
    Ri = rots_gt(:,:,i);
    
    % first self common-lines
    U1 = Ri.'*g*Ri;
    c1=[-U1(2,3)  U1(1,3)]';
    idx1 = clAngles2Ind(c1,n_theta);
    sclmatrix_gt(1,i) = idx1;
    % third self common-lines
    U2 = Ri.'*g^(n_symm-1)*Ri;
    c2=[-U2(2,3)  U2(1,3)]';
    idx2 = clAngles2Ind(c2,n_theta);
    sclmatrix_gt(2,i) = idx2;
    
    %npf_k = nonTop_npf(:,:,i);
    %correlations_selfcl(i) = npf_k(:,selfCL_matrix(1,i)).'*npf_k(:,selfCL_matrix(3,i));
end

% if strcmp(params.SCL,'GT') && params.confuse_scl
%     heads = find(round(rand(1,nImages)));
%     sclmatrix_gt(:,heads) = flipud(sclmatrix_gt(:,heads));
% end

end